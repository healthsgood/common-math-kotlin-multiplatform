package lemmingapex.trilateration

import PERFORMANCE_TEST_COUNT
import com.lemmingapex.trilateration.TrilaterationTest
import util.FastMath.sqrt
import kotlin.system.getTimeMillis

/**
 * Test class which is initialized with different predefined test cases.
 * All test cases were defined by @author scott
 *
 * @author burfi
 */
fun main() {
    val timeMillis = getTimeMillis()

    repeat(PERFORMANCE_TEST_COUNT) {
        trilateration1DExact1()
        trilateration1DExact2()
        trilateration1DInexact()
        trilateration2DExact1()
        trilateration2DExact2()
        trilateration2DExact3()
        trilateration2DExact4()
        trilateration2DExact5()
        trilateration2DInexact1()
        trilateration2DInexact2()
        trilateration2DNonIntersecting()
        trilateration2DOverIntersecting()
        trilateration3DInexact()
        trilateration4DInexact()
    }

    println("hhg-end")
    println("耗时: ${getTimeMillis() - timeMillis}")
}

fun trilateration1DExact1() {
    val positions = arrayOf(doubleArrayOf(1.0), doubleArrayOf(2.0), doubleArrayOf(3.0))
    val distances = doubleArrayOf(1.1, 0.1, 0.9)
    val expectedPosition = doubleArrayOf(2.1)
    val acceptedDelta = 0.0001
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}

fun trilateration1DExact2() {
    val positions = arrayOf(doubleArrayOf(1000.0), doubleArrayOf(2000.0), doubleArrayOf(3000.0))
    val distances = doubleArrayOf(1100.0, 100.0, 900.0)
    val expectedPosition = doubleArrayOf(2100.0)
    val acceptedDelta = 0.0001
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}

fun trilateration1DInexact() {
    val positions = arrayOf(doubleArrayOf(1000.0), doubleArrayOf(2000.0), doubleArrayOf(3000.0))
    val distances = doubleArrayOf(1110.0, 110.0, 910.0)
    val expectedPosition = doubleArrayOf(2100.0)
    val acceptedDelta = 30.0
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DExact1() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(3.0, 1.0), doubleArrayOf(2.0, 2.0))
    val distances = doubleArrayOf(1.0, 1.0, 1.0)
    val expectedPosition = doubleArrayOf(2.0, 1.0)
    val acceptedDelta = 0.0001
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DZeroDistance() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(2.0, 1.0))
    val distances = doubleArrayOf(0.0, 1.0)
    val expectedPosition = doubleArrayOf(1.0, 1.0)
    val acceptedDelta = 0.0001
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DExact2() {
    val positions = arrayOf(doubleArrayOf(0.0, 0.0), doubleArrayOf(-1.0, 0.0), doubleArrayOf(0.0, -1.0))
    val distances = doubleArrayOf(sqrt(2.0), 1.0, 1.0)
    val expectedPosition = doubleArrayOf(-1.0, -1.0)
    val acceptedDelta = 0.0001
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DExact3() {
    val positions = arrayOf(doubleArrayOf(0.0, 0.0), doubleArrayOf(1000.0, 0.0), doubleArrayOf(0.0, 1000.0))
    val distances = doubleArrayOf(sqrt(2.0) * 1000.0, 1000.0, 1000.0)
    val expectedPosition = doubleArrayOf(1000.0, 1000.0)
    val acceptedDelta = 0.0001
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DExact4() {
    val positions =
        arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(1.0, 3.0), doubleArrayOf(8.0, 8.0), doubleArrayOf(2.0, 2.0))
    val distances = doubleArrayOf(5.0, 5.0, 6.36, 3.9)
    val expectedPosition = doubleArrayOf(5.9, 2.0)
    val acceptedDelta = 0.01
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DExact5() {
    val positions = arrayOf(doubleArrayOf(5.0, -6.0), doubleArrayOf(13.0, -15.0), doubleArrayOf(21.0, -3.0))
    val distances = doubleArrayOf(8.06, 13.97, 23.32)
    val expectedPosition = doubleArrayOf(-0.6, -11.8)
    val acceptedDelta = 0.01
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DInexact1() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(3.0, 1.0), doubleArrayOf(2.0, 2.0))
    val distances = doubleArrayOf(0.9, 1.0, 1.0)
    val expectedPosition = doubleArrayOf(2.0, 1.0)
    val acceptedDelta = 0.1
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DInexact2() {
    val positions = arrayOf(
        doubleArrayOf(5.0, -6.0),
        doubleArrayOf(13.0, -15.0),
        doubleArrayOf(21.0, -3.0),
        doubleArrayOf(12.42, -21.2)
    )
    val distances = doubleArrayOf(8.06, 13.97, 23.32, 15.31)
    val expectedPosition = doubleArrayOf(-0.6, -11.8)
    val acceptedDelta = 1.0
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DNonIntersecting() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(3.0, 1.0), doubleArrayOf(2.0, 2.0))
    val distances = doubleArrayOf(0.5, 0.5, 0.5)
    val expectedPosition = doubleArrayOf(2.0, 1.0)
    val acceptedDelta = 0.25
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DOverIntersecting() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(3.0, 1.0), doubleArrayOf(2.0, 2.0))
    val distances = doubleArrayOf(2.0, 2.0, 2.0)
    val expectedPosition = doubleArrayOf(2.0, 1.0)
    val acceptedDelta = 2.0
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DDegenerateCase1() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(1.0, 1.0), doubleArrayOf(3.0, 1.0))
    val distances = doubleArrayOf(1.0, 1.0, 1.0)
    val expectedPosition = doubleArrayOf(2.0, 1.0)
    val acceptedDelta = 0.5
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DDegenerateCase2() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(1.0, 1.0), doubleArrayOf(1.0, 1.0))
    val distances = doubleArrayOf(1.0, 1.0, 1.0)
    val expectedPosition = doubleArrayOf(1.0, 1.0)
    val acceptedDelta = 0.5
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration2DUnderdertermined() {
    val positions = arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(3.0, 1.0))
    val distances = doubleArrayOf(1.0, 1.0)
    val expectedPosition = doubleArrayOf(2.0, 1.0)
    val acceptedDelta = 0.5
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration3DExact() {
    val positions =
        arrayOf(doubleArrayOf(1.0, 1.0, 1.0), doubleArrayOf(3.0, 1.0, 1.0), doubleArrayOf(2.0, 2.0, 1.0))
    val distances = doubleArrayOf(1.0, 1.0, 1.0)
    val expectedPosition = doubleArrayOf(2.0, 1.0, 1.0)
    val acceptedDelta = 0.0001
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration3DInexact() {
    val positions = arrayOf(
        doubleArrayOf(0.0, 0.0, 0.0),
        doubleArrayOf(8.84, 4.57, 12.59),
        doubleArrayOf(0.0, -8.84, 8.84),
        doubleArrayOf(10.72, -8.96, 8.84)
    )
    val distances = doubleArrayOf(8.84, 8.84, 8.84, 8.84)
    val expectedPosition = doubleArrayOf(5.2, -1.2, 7.7)
    val acceptedDelta = 1.0
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}


fun trilateration4DInexact() {
    val positions = arrayOf(
        doubleArrayOf(0.0, 0.0, 0.0, 0.0),
        doubleArrayOf(8.84, 4.57, 12.59, 9.2),
        doubleArrayOf(0.0, -8.84, 8.84, 9.2),
        doubleArrayOf(10.72, -8.96, 8.84, 9.2)
    )
    val distances = doubleArrayOf(8.84, 8.84, 8.84, 8.84)
    val expectedPosition = doubleArrayOf(5.2, -1.5, 7.7, 5.9)
    val acceptedDelta = 1.0
    TrilaterationTest(positions, distances, expectedPosition, acceptedDelta)
}
