package trilateration

import leastsquares.LeastSquaresOptimizer.Optimum
import leastsquares.LevenbergMarquardtOptimizer
import linear.RealVector

/**
 * Test class which is initialized with different predefined test cases.
 * Test was refactored from @author scott
 *
 * @author burfi
 */
class TrilaterationTest(
    private val positions: Array<DoubleArray>,
    private val distances: DoubleArray,
    private val expectedPosition: DoubleArray,
    private val acceptedDelta: Double
) {
    private var output: StringBuilder? = null
    private var linearCalculatedPosition: RealVector? = null
    private var nonLinearOptimum: Optimum? = null
    private fun testCase() {
        val trilaterationFunction = TrilaterationFunction(
            positions, distances
        )
        val lSolver = LinearLeastSquaresSolver(trilaterationFunction)
        val nlSolver = NonLinearLeastSquaresSolver(trilaterationFunction, LevenbergMarquardtOptimizer())
        linearCalculatedPosition = lSolver.solve()
        nonLinearOptimum = nlSolver.solve()
    }

    private fun outputResult() {
        output = StringBuilder()
        expectedPosition
        linearCalculatedPosition!!.toArray()
        nonLinearOptimum!!.getPoint()!!.toArray()
        output!!.append("numberOfIterations: ").append(nonLinearOptimum!!.getIterations()).append("\n")
        output!!.append("numberOfEvaluations: ").append(nonLinearOptimum!!.getEvaluations()).append("\n")
        val standardDeviation = nonLinearOptimum!!.getSigma(0.0)
        standardDeviation!!.toArray()
        output!!.append("Norm of deviation: ").append(standardDeviation.getNorm()).append("\n")
        val covarianceMatrix = nonLinearOptimum!!.getCovariances(0.0)
        output!!.append("covarianceMatrix: ").append(covarianceMatrix).append("\n")
        output.toString()
    }

    private fun compareExpectedAndCalculatedResults() {
        val calculatedPosition = nonLinearOptimum!!.getPoint()!!.toArray()
        for (i in calculatedPosition.indices) {
            expectedPosition[i]
            calculatedPosition[i]
            acceptedDelta
        }
    }

    private fun printDoubleArray(tag: String, values: DoubleArray) {
        output!!.append(tag)
        for (p in values) {
            output!!.append(p).append(" ")
        }
        output!!.append("\n")
    }

    init {
        testCase()
        outputResult()
        compareExpectedAndCalculatedResults()
    }
}