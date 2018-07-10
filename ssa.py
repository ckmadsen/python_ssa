import sys
import random
import math
from libsbml import *

def SSA(model, timeLimit, filename):
    # open file
    file = open(filename + ".csv", "w")
    header = "\"time\""
    line = "0"

	# load values from model
    values = {}
    for species in model.getListOfSpecies():
        values[species.getId()] = species.getInitialAmount()
        header += ",\"" + species.getId() + "\""
        line += "," + str(species.getInitialAmount())
    for parameter in model.getListOfParameters():
        values[parameter.getId()] = parameter.getValue()
    for reaction in model.getListOfReactions():
        for parameter in reaction.getKineticLaw().getListOfParameters():
            values[reaction.getId() + "_" + parameter.getId()] = parameter.getValue()

    # write file header and initial values
    file.write(header + "\n")
    file.write(line + "\n")

    # initialize start time
    time = 0

    # loop until time limit is reached
    while time < timeLimit:

        # initialize propensities
        propensities = {}
        for reaction in model.getListOfReactions():
            propensity = evaluateAST(reaction.getKineticLaw().getMath(), values, reaction.getId())
            for reactant in reaction.getListOfReactants():
                propensity = propensity * (values[reactant.getSpecies()] ** reactant.getStoichiometry())
            propensities[reaction.getId()] = propensity

        # evaluate the sum of all propensities
        propSum = 0
        for reaction in propensities:
            propSum += propensities[reaction]

        # draw two random numbers
        r1 = random.uniform(0, 1)
        r2 = random.uniform(0, 1)

        # determine next fire time
        fireTime = (1 / propSum) * math.log(1 / r1)

        # determine next reaction to fire
        nextReactionValue = propSum * r2
        nextReaction = None
        reactionSum = 0
        for reaction in propensities:
            reactionSum += propensities[reaction]
            if reactionSum > nextReactionValue:
                nextReaction = reaction
                break

        # update the state by firing the reaction
        time += fireTime
        if time > timeLimit:
            time = timeLimit
        for product in model.getReaction(nextReaction).getListOfProducts():
            values[product.getSpecies()] += product.getStoichiometry()
        for reactant in model.getReaction(nextReaction).getListOfReactants():
            values[reactant.getSpecies()] -= reactant.getStoichiometry()
            if values[reactant.getSpecies()] < 0:
                values[reactant.getSpecies()] = 0

        # record current state
        line = str(time)
        for species in model.getListOfSpecies():
            line += "," + str(values[species.getId()])
        file.write(line + "\n")

    # close file
    file.close()

def evaluateAST(ast, values, reactionId):
    if ast.getType() == AST_PLUS:
        return (evaluateAST(ast.getLeftChild(), values, reactionId) + evaluateAST(ast.getRightChild(), values, reactionId))
    elif ast.getType() == AST_MINUS:
        return (evaluateAST(ast.getLeftChild(), values, reactionId) - evaluateAST(ast.getRightChild(), values, reactionId))
    elif ast.getType() == AST_TIMES:
        return (evaluateAST(ast.getLeftChild(), values, reactionId) * evaluateAST(ast.getRightChild(), values, reactionId))
    elif ast.getType() == AST_DIVIDE:
        return (evaluateAST(ast.getLeftChild(), values, reactionId) / evaluateAST(ast.getRightChild(), values, reactionId))
    elif ast.getType() == AST_POWER:
        return (evaluateAST(ast.getLeftChild(), values, reactionId) ** evaluateAST(ast.getRightChild(), values, reactionId))
    elif ast.getType() == AST_FUNCTION_POWER:
        return (evaluateAST(ast.getLeftChild(), values, reactionId) ** evaluateAST(ast.getRightChild(), values, reactionId))
    elif ast.getType() == AST_INTEGER:
        return ast.getInteger()
    elif ast.getType() == AST_REAL:
        return ast.getReal()
    elif ast.getType() == AST_RATIONAL:
        return (ast.getNumerator() / ast.getDenominator())
    elif ast.getType() == AST_REAL_E:
        return (ast.getMantissa() ** ast.getExponent())
    elif ast.getType() == AST_NAME:
        if ast.getName() in values:
            return values[ast.getName()]
        else:
            return values[reactionId + "_" + ast.getName()]

def main (args):
    filename = args[1]
    timeLimit = float(args[2])
    numberOfRuns = int(args[3])
    document = readSBML(filename)
    model = document.getModel()
    for i in range(numberOfRuns):
        SSA(model, timeLimit, "run-" + str(i+1))

if __name__ == '__main__':
    main(sys.argv)