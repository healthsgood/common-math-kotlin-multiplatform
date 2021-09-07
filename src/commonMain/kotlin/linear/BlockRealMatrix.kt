package linear

import crossjvm.ArrayUtil
import exception.*
import exception.util.LocalizedFormats
import linear.exception.MatrixDimensionMismatchException
import linear.visitor.RealMatrixChangingVisitor
import linear.visitor.RealMatrixPreservingVisitor
import util.FastMath.abs
import util.FastMath.max
import util.FastMath.min
import util.FastMath.sqrt
import util.MathUtils.checkNotNull

/**
 * Cache-friendly implementation of RealMatrix using a flat arrays to store
 * square blocks of the matrix.
 *
 *
 * This implementation is specially designed to be cache-friendly. Square blocks are
 * stored as small arrays and allow efficient traversal of data both in row major direction
 * and columns major direction, one block at a time. This greatly increases performances
 * for algorithms that use crossed directions loops like multiplication or transposition.
 *
 *
 *
 * The size of square blocks is a static parameter. It may be tuned according to the cache
 * size of the target computer processor. As a rule of thumbs, it should be the largest
 * value that allows three blocks to be simultaneously cached (this is necessary for example
 * for matrix multiplication). The default value is to use 52x52 blocks which is well suited
 * for processors with 64k L1 cache (one block holds 2704 values or 21632 bytes). This value
 * could be lowered to 36x36 for processors with 32k L1 cache.
 *
 *
 *
 * The regular blocks represent [.BLOCK_SIZE] x [.BLOCK_SIZE] squares. Blocks
 * at right hand side and bottom side which may be smaller to fit matrix dimensions. The square
 * blocks are flattened in row major order in single dimension arrays which are therefore
 * [.BLOCK_SIZE]<sup>2</sup> elements long for regular blocks. The blocks are themselves
 * organized in row major order.
 *
 *
 *
 * As an example, for a block size of 52x52, a 100x60 matrix would be stored in 4 blocks.
 * Block 0 would be a double[2704] array holding the upper left 52x52 square, block 1 would be
 * a double[416] array holding the upper right 52x8 rectangle, block 2 would be a double[2496]
 * array holding the lower left 48x52 rectangle and block 3 would be a double[384] array
 * holding the lower right 48x8 rectangle.
 *
 *
 *
 * The layout complexity overhead versus simple mapping of matrices to java
 * arrays is negligible for small matrices (about 1%). The gain from cache efficiency leads
 * to up to 3-fold improvements for matrices of moderate to large size.
 *
 *
 * @since 2.0
 */
class BlockRealMatrix : AbstractRealMatrix {
    /**
     * Blocks of matrix entries.
     */
    private val blocks: Array<DoubleArray?>

    /**
     * Number of rows of the matrix.
     */
    private val rows: Int

    /**
     * Number of columns of the matrix.
     */
    private val columns: Int

    /**
     * Number of block rows of the matrix.
     */
    private val blockRows: Int

    /**
     * Number of block columns of the matrix.
     */
    private val blockColumns: Int

    /**
     * Create a new matrix with the supplied row and column dimensions.
     *
     * @param rows    the number of rows in the new matrix
     * @param columns the number of columns in the new matrix
     * @throws NotStrictlyPositiveException if row or column dimension is not
     * positive.
     */
    constructor(rows: Int, columns: Int) : super(rows, columns) {
        this.rows = rows
        this.columns = columns

        // number of blocks
        blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE
        blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE

        // allocate storage blocks, taking care of smaller ones at right and bottom
        blocks = createBlocksLayout(rows, columns)
    }

    /**
     * Create a new dense matrix copying entries from raw layout data.
     *
     * The input array *must* already be in raw layout.
     *
     * Calling this constructor is equivalent to call:
     * <pre>matrix = new BlockRealMatrix(rawData.length, rawData[0].length,
     * toBlocksLayout(rawData), false);</pre>
     *
     *
     * @param rawData data for new matrix, in raw layout
     * @throws DimensionMismatchException   if the shape of `blockData` is
     * inconsistent with block layout.
     * @throws NotStrictlyPositiveException if row or column dimension is not
     * positive.
     * @see .BlockRealMatrix
     */
    constructor(rawData: Array<DoubleArray>) : this(rawData.size, rawData[0].size, toBlocksLayout(rawData), false)

    /**
     * Create a new dense matrix copying entries from block layout data.
     *
     * The input array *must* already be in blocks layout.
     *
     * @param rows      Number of rows in the new matrix.
     * @param columns   Number of columns in the new matrix.
     * @param blockData data for new matrix
     * @param copyArray Whether the input array will be copied or referenced.
     * @throws DimensionMismatchException   if the shape of `blockData` is
     * inconsistent with block layout.
     * @throws NotStrictlyPositiveException if row or column dimension is not
     * positive.
     * @see .createBlocksLayout
     * @see .toBlocksLayout
     * @see .BlockRealMatrix
     */
    constructor(
        rows: Int, columns: Int,
        blockData: Array<DoubleArray?>, copyArray: Boolean
    ) : super(rows, columns) {
        this.rows = rows
        this.columns = columns

        // number of blocks
        blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE
        blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE
        blocks = if (copyArray) {
            // allocate storage blocks, taking care of smaller ones at right and bottom
            arrayOfNulls(blockRows * blockColumns)
        } else {
            // reference existing array
            blockData
        }
        var index = 0
        for (iBlock in 0 until blockRows) {
            val iHeight = blockHeight(iBlock)
            var jBlock = 0
            while (jBlock < blockColumns) {
                if (blockData[index]!!.size != iHeight * blockWidth(jBlock)) {
                    throw DimensionMismatchException(
                        blockData[index]!!.size,
                        iHeight * blockWidth(jBlock)
                    )
                }
                if (copyArray) {
                    blocks[index] = blockData[index]!!.copyOf()
                }
                ++jBlock
                ++index
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(NotStrictlyPositiveException::class)
    override fun createMatrix(
        rowDimension: Int,
        columnDimension: Int
    ): BlockRealMatrix? {
        return BlockRealMatrix(rowDimension, columnDimension)
    }

    /**
     * {@inheritDoc}
     */
    override fun copy(): BlockRealMatrix? {
        // create an empty matrix
        val copied = BlockRealMatrix(rows, columns)

        // copy the blocks
        for (i in blocks.indices) {
            ArrayUtil.arraycopy(blocks[i], 0, copied.blocks[i], 0, blocks[i]!!.size)
        }
        return copied
    }

    /**
     * {@inheritDoc}
     */
    @Throws(MatrixDimensionMismatchException::class)
    override fun add(m: RealMatrix?): BlockRealMatrix? {
        return try {
            add(m as BlockRealMatrix?)
        } catch (cce: ClassCastException) {
            // safety check
            MatrixUtils.checkAdditionCompatible(this, m as RealMatrix)
            val out = BlockRealMatrix(rows, columns)

            // perform addition block-wise, to ensure good cache behavior
            var blockIndex = 0
            var iBlock = 0
            while (iBlock < out.blockRows) {
                var jBlock = 0
                while (jBlock < out.blockColumns) {


                    // perform addition on the current block
                    val outBlock = out.blocks[blockIndex]
                    val tBlock = blocks[blockIndex]
                    val pStart = iBlock * BLOCK_SIZE
                    val pEnd = min(pStart + BLOCK_SIZE, rows)
                    val qStart = jBlock * BLOCK_SIZE
                    val qEnd = min(qStart + BLOCK_SIZE, columns)
                    var k = 0
                    var p = pStart
                    while (p < pEnd) {
                        var q = qStart
                        while (q < qEnd) {
                            outBlock!![k] = tBlock!![k] + m.getEntry(p, q)
                            ++k
                            ++q
                        }
                        ++p
                    }
                    // go to next block
                    ++blockIndex
                    ++jBlock
                }
                ++iBlock
            }
            out
        }
    }

    /**
     * Compute the sum of this matrix and `m`.
     *
     * @param m Matrix to be added.
     * @return `this` + m.
     * @throws MatrixDimensionMismatchException if `m` is not the same
     * size as this matrix.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun add(m: BlockRealMatrix?): BlockRealMatrix {
        // safety check
        MatrixUtils.checkAdditionCompatible(this, m as AnyMatrix)
        val out = BlockRealMatrix(rows, columns)

        // perform addition block-wise, to ensure good cache behavior
        for (blockIndex in out.blocks.indices) {
            val outBlock = out.blocks[blockIndex]
            val tBlock = blocks[blockIndex]
            val mBlock = m.blocks[blockIndex]
            for (k in outBlock!!.indices) {
                outBlock[k] = tBlock!![k] + mBlock!![k]
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(MatrixDimensionMismatchException::class)
    override fun subtract(m: RealMatrix?): BlockRealMatrix? {
        return try {
            subtract(m as BlockRealMatrix?)
        } catch (cce: ClassCastException) {
            // safety check
            MatrixUtils.checkSubtractionCompatible(this, m as AnyMatrix)
            val out = BlockRealMatrix(rows, columns)

            // perform subtraction block-wise, to ensure good cache behavior
            var blockIndex = 0
            var iBlock = 0
            while (iBlock < out.blockRows) {
                var jBlock = 0
                while (jBlock < out.blockColumns) {


                    // perform subtraction on the current block
                    val outBlock = out.blocks[blockIndex]
                    val tBlock = blocks[blockIndex]
                    val pStart = iBlock * BLOCK_SIZE
                    val pEnd = min(pStart + BLOCK_SIZE, rows)
                    val qStart = jBlock * BLOCK_SIZE
                    val qEnd = min(qStart + BLOCK_SIZE, columns)
                    var k = 0
                    var p = pStart
                    while (p < pEnd) {
                        var q = qStart
                        while (q < qEnd) {
                            outBlock!![k] = tBlock!![k] - m.getEntry(p, q)
                            ++k
                            ++q
                        }
                        ++p
                    }
                    // go to next block
                    ++blockIndex
                    ++jBlock
                }
                ++iBlock
            }
            out
        }
    }

    /**
     * Subtract `m` from this matrix.
     *
     * @param m Matrix to be subtracted.
     * @return `this` - m.
     * @throws MatrixDimensionMismatchException if `m` is not the
     * same size as this matrix.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun subtract(m: BlockRealMatrix?): BlockRealMatrix {
        // safety check
        MatrixUtils.checkSubtractionCompatible(this, m as AnyMatrix)
        val out = BlockRealMatrix(rows, columns)

        // perform subtraction block-wise, to ensure good cache behavior
        for (blockIndex in out.blocks.indices) {
            val outBlock = out.blocks[blockIndex]
            val tBlock = blocks[blockIndex]
            val mBlock = m.blocks[blockIndex]
            for (k in outBlock!!.indices) {
                outBlock[k] = tBlock!![k] - mBlock!![k]
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun scalarAdd(d: Double): BlockRealMatrix? {
        val out = BlockRealMatrix(rows, columns)

        // perform subtraction block-wise, to ensure good cache behavior
        for (blockIndex in out.blocks.indices) {
            val outBlock = out.blocks[blockIndex]
            val tBlock = blocks[blockIndex]
            for (k in outBlock!!.indices) {
                outBlock[k] = tBlock!![k] + d
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun scalarMultiply(d: Double): RealMatrix? {
        val out = BlockRealMatrix(rows, columns)

        // perform subtraction block-wise, to ensure good cache behavior
        for (blockIndex in out.blocks.indices) {
            val outBlock = out.blocks[blockIndex]
            val tBlock = blocks[blockIndex]
            for (k in outBlock!!.indices) {
                outBlock[k] = tBlock!![k] * d
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun multiply(m: RealMatrix?): BlockRealMatrix? {
        return try {
            multiply(m as BlockRealMatrix?)
        } catch (cce: ClassCastException) {
            // safety check
            MatrixUtils.checkMultiplicationCompatible(this, m as AnyMatrix)
            val out = BlockRealMatrix(rows, m.getColumnDimension())

            // perform multiplication block-wise, to ensure good cache behavior
            var blockIndex = 0
            var iBlock = 0
            while (iBlock < out.blockRows) {
                val pStart = iBlock * BLOCK_SIZE
                val pEnd = min(pStart + BLOCK_SIZE, rows)
                var jBlock = 0
                while (jBlock < out.blockColumns) {
                    val qStart = jBlock * BLOCK_SIZE
                    val qEnd = min(qStart + BLOCK_SIZE, m.getColumnDimension())

                    // select current block
                    val outBlock = out.blocks[blockIndex]

                    // perform multiplication on current block
                    var kBlock = 0
                    while (kBlock < blockColumns) {
                        val kWidth = blockWidth(kBlock)
                        val tBlock = blocks[iBlock * blockColumns + kBlock]
                        val rStart = kBlock * BLOCK_SIZE
                        var k = 0
                        var p = pStart
                        while (p < pEnd) {
                            val lStart = (p - pStart) * kWidth
                            val lEnd = lStart + kWidth
                            var q = qStart
                            while (q < qEnd) {
                                var sum = 0.0
                                var r = rStart
                                var l = lStart
                                while (l < lEnd) {
                                    sum += tBlock!![l] * m.getEntry(r, q)
                                    ++r
                                    ++l
                                }
                                outBlock!![k] += sum
                                ++k
                                ++q
                            }
                            ++p
                        }
                        ++kBlock
                    }
                    // go to next block
                    ++blockIndex
                    ++jBlock
                }
                ++iBlock
            }
            out
        }
    }

    /**
     * Returns the result of postmultiplying this by `m`.
     *
     * @param m Matrix to postmultiply by.
     * @return `this` * m.
     * @throws DimensionMismatchException if the matrices are not compatible.
     */
    @Throws(DimensionMismatchException::class)
    fun multiply(m: BlockRealMatrix?): BlockRealMatrix {
        // safety check
        MatrixUtils.checkMultiplicationCompatible(this, m as AnyMatrix)
        val out = BlockRealMatrix(rows, m.columns)

        // perform multiplication block-wise, to ensure good cache behavior
        var blockIndex = 0
        for (iBlock in 0 until out.blockRows) {
            val pStart = iBlock * BLOCK_SIZE
            val pEnd = min(pStart + BLOCK_SIZE, rows)
            for (jBlock in 0 until out.blockColumns) {
                val jWidth = out.blockWidth(jBlock)
                val jWidth2 = jWidth + jWidth
                val jWidth3 = jWidth2 + jWidth
                val jWidth4 = jWidth3 + jWidth

                // select current block
                val outBlock = out.blocks[blockIndex]

                // perform multiplication on current block
                for (kBlock in 0 until blockColumns) {
                    val kWidth = blockWidth(kBlock)
                    val tBlock = blocks[iBlock * blockColumns + kBlock]
                    val mBlock = m.blocks[kBlock * m.blockColumns + jBlock]
                    var k = 0
                    for (p in pStart until pEnd) {
                        val lStart = (p - pStart) * kWidth
                        val lEnd = lStart + kWidth
                        for (nStart in 0 until jWidth) {
                            var sum = 0.0
                            var l = lStart
                            var n = nStart
                            while (l < lEnd - 3) {
                                sum += tBlock!![l] * mBlock!![n] + tBlock[l + 1] * mBlock[n + jWidth] + tBlock[l + 2] * mBlock[n + jWidth2] + tBlock[l + 3] * mBlock[n + jWidth3]
                                l += 4
                                n += jWidth4
                            }
                            while (l < lEnd) {
                                sum += tBlock!![l++] * mBlock!![n]
                                n += jWidth
                            }
                            outBlock!![k] += sum
                            ++k
                        }
                    }
                }
                // go to next block
                ++blockIndex
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun getData(): Array<DoubleArray?>? {
        val data = Array<DoubleArray?>(getRowDimension()) { DoubleArray(getColumnDimension()) }
        val lastColumns = columns - (blockColumns - 1) * BLOCK_SIZE
        for (iBlock in 0 until blockRows) {
            val pStart = iBlock * BLOCK_SIZE
            val pEnd = min(pStart + BLOCK_SIZE, rows)
            var regularPos = 0
            var lastPos = 0
            for (p in pStart until pEnd) {
                val dataP = data[p]
                var blockIndex = iBlock * blockColumns
                var dataPos = 0
                for (jBlock in 0 until blockColumns - 1) {
                    ArrayUtil.arraycopy(blocks[blockIndex++], regularPos, dataP, dataPos, BLOCK_SIZE)
                    dataPos += BLOCK_SIZE
                }
                ArrayUtil.arraycopy(blocks[blockIndex], lastPos, dataP, dataPos, lastColumns)
                regularPos += BLOCK_SIZE
                lastPos += lastColumns
            }
        }
        return data
    }

    /**
     * {@inheritDoc}
     */
    override fun getNorm(): Double {
        val colSums = DoubleArray(BLOCK_SIZE)
        var maxColSum = 0.0
        for (jBlock in 0 until blockColumns) {
            val jWidth = blockWidth(jBlock)
            ArrayUtil.fill(colSums, 0, jWidth, 0.0)
            for (iBlock in 0 until blockRows) {
                val iHeight = blockHeight(iBlock)
                val block = blocks[iBlock * blockColumns + jBlock]
                for (j in 0 until jWidth) {
                    var sum = 0.0
                    for (i in 0 until iHeight) {
                        sum += abs(block!![i * jWidth + j])
                    }
                    colSums[j] += sum
                }
            }
            for (j in 0 until jWidth) {
                maxColSum = max(maxColSum, colSums[j])
            }
        }
        return maxColSum
    }

    /**
     * {@inheritDoc}
     */
    override fun getFrobeniusNorm(): Double {
        var sum2 = 0.0
        for (blockIndex in blocks.indices) {
            for (entry in blocks[blockIndex]!!) {
                sum2 += entry * entry
            }
        }
        return sqrt(sum2)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun getSubMatrix(
        startRow: Int, endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): BlockRealMatrix? {
        // safety checks
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)

        // create the output matrix
        val out = BlockRealMatrix(endRow - startRow + 1, endColumn - startColumn + 1)

        // compute blocks shifts
        val blockStartRow = startRow / BLOCK_SIZE
        val rowsShift = startRow % BLOCK_SIZE
        val blockStartColumn = startColumn / BLOCK_SIZE
        val columnsShift = startColumn % BLOCK_SIZE

        // perform extraction block-wise, to ensure good cache behavior
        var pBlock = blockStartRow
        for (iBlock in 0 until out.blockRows) {
            val iHeight = out.blockHeight(iBlock)
            var qBlock = blockStartColumn
            for (jBlock in 0 until out.blockColumns) {
                val jWidth = out.blockWidth(jBlock)

                // handle one block of the output matrix
                val outIndex = iBlock * out.blockColumns + jBlock
                val outBlock = out.blocks[outIndex]
                val index = pBlock * blockColumns + qBlock
                val width = blockWidth(qBlock)
                val heightExcess = iHeight + rowsShift - BLOCK_SIZE
                val widthExcess = jWidth + columnsShift - BLOCK_SIZE
                if (heightExcess > 0) {
                    // the submatrix block spans on two blocks rows from the original matrix
                    if (widthExcess > 0) {
                        // the submatrix block spans on two blocks columns from the original matrix
                        val width2 = blockWidth(qBlock + 1)
                        copyBlockPart(
                            blocks[index], width,
                            rowsShift, BLOCK_SIZE,
                            columnsShift, BLOCK_SIZE,
                            outBlock, jWidth, 0, 0
                        )
                        copyBlockPart(
                            blocks[index + 1], width2,
                            rowsShift, BLOCK_SIZE,
                            0, widthExcess,
                            outBlock, jWidth, 0, jWidth - widthExcess
                        )
                        copyBlockPart(
                            blocks[index + blockColumns], width,
                            0, heightExcess,
                            columnsShift, BLOCK_SIZE,
                            outBlock, jWidth, iHeight - heightExcess, 0
                        )
                        copyBlockPart(
                            blocks[index + blockColumns + 1], width2,
                            0, heightExcess,
                            0, widthExcess,
                            outBlock, jWidth, iHeight - heightExcess, jWidth - widthExcess
                        )
                    } else {
                        // the submatrix block spans on one block column from the original matrix
                        copyBlockPart(
                            blocks[index], width,
                            rowsShift, BLOCK_SIZE,
                            columnsShift, jWidth + columnsShift,
                            outBlock, jWidth, 0, 0
                        )
                        copyBlockPart(
                            blocks[index + blockColumns], width,
                            0, heightExcess,
                            columnsShift, jWidth + columnsShift,
                            outBlock, jWidth, iHeight - heightExcess, 0
                        )
                    }
                } else {
                    // the submatrix block spans on one block row from the original matrix
                    if (widthExcess > 0) {
                        // the submatrix block spans on two blocks columns from the original matrix
                        val width2 = blockWidth(qBlock + 1)
                        copyBlockPart(
                            blocks[index], width,
                            rowsShift, iHeight + rowsShift,
                            columnsShift, BLOCK_SIZE,
                            outBlock, jWidth, 0, 0
                        )
                        copyBlockPart(
                            blocks[index + 1], width2,
                            rowsShift, iHeight + rowsShift,
                            0, widthExcess,
                            outBlock, jWidth, 0, jWidth - widthExcess
                        )
                    } else {
                        // the submatrix block spans on one block column from the original matrix
                        copyBlockPart(
                            blocks[index], width,
                            rowsShift, iHeight + rowsShift,
                            columnsShift, jWidth + columnsShift,
                            outBlock, jWidth, 0, 0
                        )
                    }
                }
                ++qBlock
            }
            ++pBlock
        }
        return out
    }

    /**
     * Copy a part of a block into another one
     *
     * This method can be called only when the specified part fits in both
     * blocks, no verification is done here.
     *
     * @param srcBlock       source block
     * @param srcWidth       source block width ([.BLOCK_SIZE] or smaller)
     * @param srcStartRow    start row in the source block
     * @param srcEndRow      end row (exclusive) in the source block
     * @param srcStartColumn start column in the source block
     * @param srcEndColumn   end column (exclusive) in the source block
     * @param dstBlock       destination block
     * @param dstWidth       destination block width ([.BLOCK_SIZE] or smaller)
     * @param dstStartRow    start row in the destination block
     * @param dstStartColumn start column in the destination block
     */
    private fun copyBlockPart(
        srcBlock: DoubleArray?, srcWidth: Int,
        srcStartRow: Int, srcEndRow: Int,
        srcStartColumn: Int, srcEndColumn: Int,
        dstBlock: DoubleArray?, dstWidth: Int,
        dstStartRow: Int, dstStartColumn: Int
    ) {
        val length = srcEndColumn - srcStartColumn
        var srcPos = srcStartRow * srcWidth + srcStartColumn
        var dstPos = dstStartRow * dstWidth + dstStartColumn
        for (srcRow in srcStartRow until srcEndRow) {
            ArrayUtil.arraycopy(srcBlock, srcPos, dstBlock, dstPos, length)
            srcPos += srcWidth
            dstPos += dstWidth
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(
        OutOfRangeException::class,
        NoDataException::class,
        NullArgumentException::class,
        DimensionMismatchException::class
    )
    override fun setSubMatrix(
        subMatrix: Array<DoubleArray?>?,
        row: Int,
        column: Int
    ) {
        // safety checks
        checkNotNull(subMatrix)
        val refLength: Int = subMatrix!![0]!!.size
        if (refLength == 0) {
            throw NoDataException(LocalizedFormats.AT_LEAST_ONE_COLUMN)
        }
        val endRow = row + subMatrix.size - 1
        val endColumn = column + refLength - 1
        MatrixUtils.checkSubMatrixIndex(this, row, endRow, column, endColumn)
        for (subRow in subMatrix) {
            if (subRow!!.size != refLength) {
                throw DimensionMismatchException(refLength, subRow.size)
            }
        }

        // compute blocks bounds
        val blockStartRow = row / BLOCK_SIZE
        val blockEndRow = (endRow + BLOCK_SIZE) / BLOCK_SIZE
        val blockStartColumn = column / BLOCK_SIZE
        val blockEndColumn = (endColumn + BLOCK_SIZE) / BLOCK_SIZE

        // perform copy block-wise, to ensure good cache behavior
        for (iBlock in blockStartRow until blockEndRow) {
            val iHeight = blockHeight(iBlock)
            val firstRow = iBlock * BLOCK_SIZE
            val iStart = max(row, firstRow)
            val iEnd = min(endRow + 1, firstRow + iHeight)
            for (jBlock in blockStartColumn until blockEndColumn) {
                val jWidth = blockWidth(jBlock)
                val firstColumn = jBlock * BLOCK_SIZE
                val jStart = max(column, firstColumn)
                val jEnd = min(endColumn + 1, firstColumn + jWidth)
                val jLength = jEnd - jStart

                // handle one block, row by row
                val block = blocks[iBlock * blockColumns + jBlock]
                for (i in iStart until iEnd) {
                    ArrayUtil.arraycopy(
                        subMatrix[i - row], jStart - column,
                        block, (i - firstRow) * jWidth + (jStart - firstColumn),
                        jLength
                    )
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getRowMatrix(row: Int): BlockRealMatrix? {
        MatrixUtils.checkRowIndex(this, row)
        val out = BlockRealMatrix(1, columns)

        // perform copy block-wise, to ensure good cache behavior
        val iBlock = row / BLOCK_SIZE
        val iRow = row - iBlock * BLOCK_SIZE
        var outBlockIndex = 0
        var outIndex = 0
        var outBlock = out.blocks[outBlockIndex]
        for (jBlock in 0 until blockColumns) {
            val jWidth = blockWidth(jBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            val available = outBlock!!.size - outIndex
            if (jWidth > available) {
                ArrayUtil.arraycopy(block, iRow * jWidth, outBlock, outIndex, available)
                outBlock = out.blocks[++outBlockIndex]
                ArrayUtil.arraycopy(block, iRow * jWidth, outBlock, 0, jWidth - available)
                outIndex = jWidth - available
            } else {
                ArrayUtil.arraycopy(block, iRow * jWidth, outBlock, outIndex, jWidth)
                outIndex += jWidth
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setRowMatrix(row: Int, matrix: RealMatrix?) {
        try {
            setRowMatrix(row, matrix as BlockRealMatrix?)
        } catch (cce: ClassCastException) {
            super.setRowMatrix(row, matrix)
        }
    }

    /**
     * Sets the entries in row number `row`
     * as a row matrix.  Row indices start at 0.
     *
     * @param row    the row to be set
     * @param matrix row matrix (must have one row and the same number of columns
     * as the instance)
     * @throws OutOfRangeException              if the specified row index is invalid.
     * @throws MatrixDimensionMismatchException if the matrix dimensions do
     * not match one instance row.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setRowMatrix(row: Int, matrix: BlockRealMatrix?) {
        MatrixUtils.checkRowIndex(this, row)
        val nCols = getColumnDimension()
        if (matrix!!.getRowDimension() != 1 ||
            matrix.getColumnDimension() != nCols
        ) {
            throw MatrixDimensionMismatchException(
                matrix.getRowDimension(),
                matrix.getColumnDimension(),
                1, nCols
            )
        }

        // perform copy block-wise, to ensure good cache behavior
        val iBlock = row / BLOCK_SIZE
        val iRow = row - iBlock * BLOCK_SIZE
        var mBlockIndex = 0
        var mIndex = 0
        var mBlock = matrix.blocks[mBlockIndex]
        for (jBlock in 0 until blockColumns) {
            val jWidth = blockWidth(jBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            val available = mBlock!!.size - mIndex
            if (jWidth > available) {
                ArrayUtil.arraycopy(mBlock, mIndex, block, iRow * jWidth, available)
                mBlock = matrix.blocks[++mBlockIndex]
                ArrayUtil.arraycopy(mBlock, 0, block, iRow * jWidth, jWidth - available)
                mIndex = jWidth - available
            } else {
                ArrayUtil.arraycopy(mBlock, mIndex, block, iRow * jWidth, jWidth)
                mIndex += jWidth
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getColumnMatrix(column: Int): BlockRealMatrix? {
        MatrixUtils.checkColumnIndex(this, column)
        val out = BlockRealMatrix(rows, 1)

        // perform copy block-wise, to ensure good cache behavior
        val jBlock = column / BLOCK_SIZE
        val jColumn = column - jBlock * BLOCK_SIZE
        val jWidth = blockWidth(jBlock)
        var outBlockIndex = 0
        var outIndex = 0
        var outBlock = out.blocks[outBlockIndex]
        for (iBlock in 0 until blockRows) {
            val iHeight = blockHeight(iBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            for (i in 0 until iHeight) {
                if (outIndex >= outBlock!!.size) {
                    outBlock = out.blocks[++outBlockIndex]
                    outIndex = 0
                }
                outBlock!![outIndex++] = block!![i * jWidth + jColumn]
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setColumnMatrix(column: Int, matrix: RealMatrix?) {
        try {
            setColumnMatrix(column, matrix as BlockRealMatrix?)
        } catch (cce: ClassCastException) {
            super.setColumnMatrix(column, matrix)
        }
    }

    /**
     * Sets the entries in column number `column`
     * as a column matrix.  Column indices start at 0.
     *
     * @param column the column to be set
     * @param matrix column matrix (must have one column and the same number of rows
     * as the instance)
     * @throws OutOfRangeException              if the specified column index is invalid.
     * @throws MatrixDimensionMismatchException if the matrix dimensions do
     * not match one instance column.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setColumnMatrix(column: Int, matrix: BlockRealMatrix?) {
        MatrixUtils.checkColumnIndex(this, column)
        val nRows = getRowDimension()
        if (matrix!!.getRowDimension() != nRows ||
            matrix.getColumnDimension() != 1
        ) {
            throw MatrixDimensionMismatchException(
                matrix.getRowDimension(),
                matrix.getColumnDimension(),
                nRows, 1
            )
        }

        // perform copy block-wise, to ensure good cache behavior
        val jBlock = column / BLOCK_SIZE
        val jColumn = column - jBlock * BLOCK_SIZE
        val jWidth = blockWidth(jBlock)
        var mBlockIndex = 0
        var mIndex = 0
        var mBlock = matrix.blocks[mBlockIndex]
        for (iBlock in 0 until blockRows) {
            val iHeight = blockHeight(iBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            for (i in 0 until iHeight) {
                if (mIndex >= mBlock!!.size) {
                    mBlock = matrix.blocks[++mBlockIndex]
                    mIndex = 0
                }
                block!![i * jWidth + jColumn] = mBlock!![mIndex++]
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getRowVector(row: Int): RealVector? {
        MatrixUtils.checkRowIndex(this, row)
        val outData = DoubleArray(columns)

        // perform copy block-wise, to ensure good cache behavior
        val iBlock = row / BLOCK_SIZE
        val iRow = row - iBlock * BLOCK_SIZE
        var outIndex = 0
        for (jBlock in 0 until blockColumns) {
            val jWidth = blockWidth(jBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            ArrayUtil.arraycopy(block, iRow * jWidth, outData, outIndex, jWidth)
            outIndex += jWidth
        }
        return ArrayRealVector(outData, false)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setRowVector(row: Int, vector: RealVector?) {
        try {
            setRow(row, (vector as ArrayRealVector?)!!.dataRef)
        } catch (cce: ClassCastException) {
            super.setRowVector(row, vector)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getColumnVector(column: Int): RealVector? {
        MatrixUtils.checkColumnIndex(this, column)
        val outData = DoubleArray(rows)

        // perform copy block-wise, to ensure good cache behavior
        val jBlock = column / BLOCK_SIZE
        val jColumn = column - jBlock * BLOCK_SIZE
        val jWidth = blockWidth(jBlock)
        var outIndex = 0
        for (iBlock in 0 until blockRows) {
            val iHeight = blockHeight(iBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            for (i in 0 until iHeight) {
                outData[outIndex++] = block!![i * jWidth + jColumn]
            }
        }
        return ArrayRealVector(outData, false)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setColumnVector(column: Int, vector: RealVector?) {
        try {
            setColumn(column, (vector as ArrayRealVector?)!!.dataRef)
        } catch (cce: ClassCastException) {
            super.setColumnVector(column, vector)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getRow(row: Int): DoubleArray? {
        MatrixUtils.checkRowIndex(this, row)
        val out = DoubleArray(columns)

        // perform copy block-wise, to ensure good cache behavior
        val iBlock = row / BLOCK_SIZE
        val iRow = row - iBlock * BLOCK_SIZE
        var outIndex = 0
        for (jBlock in 0 until blockColumns) {
            val jWidth = blockWidth(jBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            ArrayUtil.arraycopy(block, iRow * jWidth, out, outIndex, jWidth)
            outIndex += jWidth
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setRow(row: Int, array: DoubleArray?) {
        MatrixUtils.checkRowIndex(this, row)
        val nCols = getColumnDimension()
        if (array!!.size != nCols) {
            throw MatrixDimensionMismatchException(1, array.size, 1, nCols)
        }

        // perform copy block-wise, to ensure good cache behavior
        val iBlock = row / BLOCK_SIZE
        val iRow = row - iBlock * BLOCK_SIZE
        var outIndex = 0
        for (jBlock in 0 until blockColumns) {
            val jWidth = blockWidth(jBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            ArrayUtil.arraycopy(array, outIndex, block, iRow * jWidth, jWidth)
            outIndex += jWidth
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getColumn(column: Int): DoubleArray? {
        MatrixUtils.checkColumnIndex(this, column)
        val out = DoubleArray(rows)

        // perform copy block-wise, to ensure good cache behavior
        val jBlock = column / BLOCK_SIZE
        val jColumn = column - jBlock * BLOCK_SIZE
        val jWidth = blockWidth(jBlock)
        var outIndex = 0
        for (iBlock in 0 until blockRows) {
            val iHeight = blockHeight(iBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            for (i in 0 until iHeight) {
                out[outIndex++] = block!![i * jWidth + jColumn]
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setColumn(column: Int, array: DoubleArray?) {
        MatrixUtils.checkColumnIndex(this, column)
        val nRows = getRowDimension()
        if (array!!.size != nRows) {
            throw MatrixDimensionMismatchException(array.size, 1, nRows, 1)
        }

        // perform copy block-wise, to ensure good cache behavior
        val jBlock = column / BLOCK_SIZE
        val jColumn = column - jBlock * BLOCK_SIZE
        val jWidth = blockWidth(jBlock)
        var outIndex = 0
        for (iBlock in 0 until blockRows) {
            val iHeight = blockHeight(iBlock)
            val block = blocks[iBlock * blockColumns + jBlock]
            for (i in 0 until iHeight) {
                block!![i * jWidth + jColumn] = array[outIndex++]
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    override fun getEntry(row: Int, column: Int): Double {
        MatrixUtils.checkMatrixIndex(this, row, column)
        val iBlock = row / BLOCK_SIZE
        val jBlock = column / BLOCK_SIZE
        val k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE)
        return blocks[iBlock * blockColumns + jBlock]!![k]
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun setEntry(row: Int, column: Int, value: Double) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        val iBlock = row / BLOCK_SIZE
        val jBlock = column / BLOCK_SIZE
        val k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE)
        blocks[iBlock * blockColumns + jBlock]!![k] = value
    }

    /**
     * {@inheritDoc}
     */
    override fun addToEntry(
        row: Int, column: Int,
        increment: Double
    ) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        val iBlock = row / BLOCK_SIZE
        val jBlock = column / BLOCK_SIZE
        val k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE)
        blocks[iBlock * blockColumns + jBlock]!![k] += increment
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun multiplyEntry(
        row: Int, column: Int,
        factor: Double
    ) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        val iBlock = row / BLOCK_SIZE
        val jBlock = column / BLOCK_SIZE
        val k = (row - iBlock * BLOCK_SIZE) * blockWidth(jBlock) +
                (column - jBlock * BLOCK_SIZE)
        blocks[iBlock * blockColumns + jBlock]!![k] *= factor
    }

    /**
     * {@inheritDoc}
     */
    override fun transpose(): BlockRealMatrix? {
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        val out = BlockRealMatrix(nCols, nRows)

        // perform transpose block-wise, to ensure good cache behavior
        var blockIndex = 0
        for (iBlock in 0 until blockColumns) {
            for (jBlock in 0 until blockRows) {
                // transpose current block
                val outBlock = out.blocks[blockIndex]
                val tBlock = blocks[jBlock * blockColumns + iBlock]
                val pStart = iBlock * BLOCK_SIZE
                val pEnd = min(pStart + BLOCK_SIZE, columns)
                val qStart = jBlock * BLOCK_SIZE
                val qEnd = min(qStart + BLOCK_SIZE, rows)
                var k = 0
                for (p in pStart until pEnd) {
                    val lInc = pEnd - pStart
                    var l = p - pStart
                    for (q in qStart until qEnd) {
                        outBlock!![k] = tBlock!![l]
                        ++k
                        l += lInc
                    }
                }
                // go to next block
                ++blockIndex
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun getRowDimension(): Int {
        return rows
    }

    /**
     * {@inheritDoc}
     */
    override fun getColumnDimension(): Int {
        return columns
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun operate(v: DoubleArray?): DoubleArray? {
        if (v!!.size != columns) {
            throw DimensionMismatchException(v.size, columns)
        }
        val out = DoubleArray(rows)

        // perform multiplication block-wise, to ensure good cache behavior
        for (iBlock in 0 until blockRows) {
            val pStart = iBlock * BLOCK_SIZE
            val pEnd = min(pStart + BLOCK_SIZE, rows)
            for (jBlock in 0 until blockColumns) {
                val block = blocks[iBlock * blockColumns + jBlock]
                val qStart = jBlock * BLOCK_SIZE
                val qEnd = min(qStart + BLOCK_SIZE, columns)
                var k = 0
                for (p in pStart until pEnd) {
                    var sum = 0.0
                    var q = qStart
                    while (q < qEnd - 3) {
                        sum += block!![k] * v[q] + block[k + 1] * v[q + 1] + block[k + 2] * v[q + 2] + block[k + 3] * v[q + 3]
                        k += 4
                        q += 4
                    }
                    while (q < qEnd) {
                        sum += block!![k++] * v[q++]
                    }
                    out[p] += sum
                }
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun preMultiply(v: DoubleArray): DoubleArray? {
        if (v.size != rows) {
            throw DimensionMismatchException(v.size, rows)
        }
        val out = DoubleArray(columns)

        // perform multiplication block-wise, to ensure good cache behavior
        for (jBlock in 0 until blockColumns) {
            val jWidth = blockWidth(jBlock)
            val jWidth2 = jWidth + jWidth
            val jWidth3 = jWidth2 + jWidth
            val jWidth4 = jWidth3 + jWidth
            val qStart = jBlock * BLOCK_SIZE
            val qEnd = min(qStart + BLOCK_SIZE, columns)
            for (iBlock in 0 until blockRows) {
                val block = blocks[iBlock * blockColumns + jBlock]
                val pStart = iBlock * BLOCK_SIZE
                val pEnd = min(pStart + BLOCK_SIZE, rows)
                for (q in qStart until qEnd) {
                    var k = q - qStart
                    var sum = 0.0
                    var p = pStart
                    while (p < pEnd - 3) {
                        sum += block!![k] * v[p] + block[k + jWidth] * v[p + 1] + block[k + jWidth2] * v[p + 2] + block[k + jWidth3] * v[p + 3]
                        k += jWidth4
                        p += 4
                    }
                    while (p < pEnd) {
                        sum += block!![k] * v[p++]
                        k += jWidth
                    }
                    out[q] += sum
                }
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInRowOrder(visitor: RealMatrixChangingVisitor): Double {
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (iBlock in 0 until blockRows) {
            val pStart = iBlock * BLOCK_SIZE
            val pEnd = min(pStart + BLOCK_SIZE, rows)
            for (p in pStart until pEnd) {
                for (jBlock in 0 until blockColumns) {
                    val jWidth = blockWidth(jBlock)
                    val qStart = jBlock * BLOCK_SIZE
                    val qEnd = min(qStart + BLOCK_SIZE, columns)
                    val block = blocks[iBlock * blockColumns + jBlock]
                    var k = (p - pStart) * jWidth
                    for (q in qStart until qEnd) {
                        block!![k] = visitor.visit(p, q, block[k])
                        ++k
                    }
                }
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInRowOrder(visitor: RealMatrixPreservingVisitor?): Double {
        visitor!!.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (iBlock in 0 until blockRows) {
            val pStart = iBlock * BLOCK_SIZE
            val pEnd = min(pStart + BLOCK_SIZE, rows)
            for (p in pStart until pEnd) {
                for (jBlock in 0 until blockColumns) {
                    val jWidth = blockWidth(jBlock)
                    val qStart = jBlock * BLOCK_SIZE
                    val qEnd = min(qStart + BLOCK_SIZE, columns)
                    val block = blocks[iBlock * blockColumns + jBlock]
                    var k = (p - pStart) * jWidth
                    for (q in qStart until qEnd) {
                        visitor.visit(p, q, block!![k])
                        ++k
                    }
                }
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInRowOrder(
        visitor: RealMatrixChangingVisitor,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(rows, columns, startRow, endRow, startColumn, endColumn)
        for (iBlock in startRow / BLOCK_SIZE until 1 + endRow / BLOCK_SIZE) {
            val p0 = iBlock * BLOCK_SIZE
            val pStart = max(startRow, p0)
            val pEnd = min((iBlock + 1) * BLOCK_SIZE, 1 + endRow)
            for (p in pStart until pEnd) {
                for (jBlock in startColumn / BLOCK_SIZE until 1 + endColumn / BLOCK_SIZE) {
                    val jWidth = blockWidth(jBlock)
                    val q0 = jBlock * BLOCK_SIZE
                    val qStart = max(startColumn, q0)
                    val qEnd = min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn)
                    val block = blocks[iBlock * blockColumns + jBlock]
                    var k = (p - p0) * jWidth + qStart - q0
                    for (q in qStart until qEnd) {
                        block!![k] = visitor.visit(p, q, block[k])
                        ++k
                    }
                }
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInRowOrder(
        visitor: RealMatrixPreservingVisitor?,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor!!.start(rows, columns, startRow, endRow, startColumn, endColumn)
        for (iBlock in startRow / BLOCK_SIZE until 1 + endRow / BLOCK_SIZE) {
            val p0 = iBlock * BLOCK_SIZE
            val pStart = max(startRow, p0)
            val pEnd = min((iBlock + 1) * BLOCK_SIZE, 1 + endRow)
            for (p in pStart until pEnd) {
                for (jBlock in startColumn / BLOCK_SIZE until 1 + endColumn / BLOCK_SIZE) {
                    val jWidth = blockWidth(jBlock)
                    val q0 = jBlock * BLOCK_SIZE
                    val qStart = max(startColumn, q0)
                    val qEnd = min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn)
                    val block = blocks[iBlock * blockColumns + jBlock]
                    var k = (p - p0) * jWidth + qStart - q0
                    for (q in qStart until qEnd) {
                        visitor.visit(p, q, block!![k])
                        ++k
                    }
                }
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInOptimizedOrder(visitor: RealMatrixChangingVisitor): Double {
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        var blockIndex = 0
        for (iBlock in 0 until blockRows) {
            val pStart = iBlock * BLOCK_SIZE
            val pEnd = min(pStart + BLOCK_SIZE, rows)
            for (jBlock in 0 until blockColumns) {
                val qStart = jBlock * BLOCK_SIZE
                val qEnd = min(qStart + BLOCK_SIZE, columns)
                val block = blocks[blockIndex]
                var k = 0
                for (p in pStart until pEnd) {
                    for (q in qStart until qEnd) {
                        block!![k] = visitor.visit(p, q, block[k])
                        ++k
                    }
                }
                ++blockIndex
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInOptimizedOrder(visitor: RealMatrixPreservingVisitor?): Double {
        visitor!!.start(rows, columns, 0, rows - 1, 0, columns - 1)
        var blockIndex = 0
        for (iBlock in 0 until blockRows) {
            val pStart = iBlock * BLOCK_SIZE
            val pEnd = min(pStart + BLOCK_SIZE, rows)
            for (jBlock in 0 until blockColumns) {
                val qStart = jBlock * BLOCK_SIZE
                val qEnd = min(qStart + BLOCK_SIZE, columns)
                val block = blocks[blockIndex]
                var k = 0
                for (p in pStart until pEnd) {
                    for (q in qStart until qEnd) {
                        visitor.visit(p, q, block!![k])
                        ++k
                    }
                }
                ++blockIndex
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInOptimizedOrder(
        visitor: RealMatrixChangingVisitor,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(rows, columns, startRow, endRow, startColumn, endColumn)
        for (iBlock in startRow / BLOCK_SIZE until 1 + endRow / BLOCK_SIZE) {
            val p0 = iBlock * BLOCK_SIZE
            val pStart = max(startRow, p0)
            val pEnd = min((iBlock + 1) * BLOCK_SIZE, 1 + endRow)
            for (jBlock in startColumn / BLOCK_SIZE until 1 + endColumn / BLOCK_SIZE) {
                val jWidth = blockWidth(jBlock)
                val q0 = jBlock * BLOCK_SIZE
                val qStart = max(startColumn, q0)
                val qEnd = min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn)
                val block = blocks[iBlock * blockColumns + jBlock]
                for (p in pStart until pEnd) {
                    var k = (p - p0) * jWidth + qStart - q0
                    for (q in qStart until qEnd) {
                        block!![k] = visitor.visit(p, q, block[k])
                        ++k
                    }
                }
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInOptimizedOrder(
        visitor: RealMatrixPreservingVisitor?,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor!!.start(rows, columns, startRow, endRow, startColumn, endColumn)
        for (iBlock in startRow / BLOCK_SIZE until 1 + endRow / BLOCK_SIZE) {
            val p0 = iBlock * BLOCK_SIZE
            val pStart = max(startRow, p0)
            val pEnd = min((iBlock + 1) * BLOCK_SIZE, 1 + endRow)
            for (jBlock in startColumn / BLOCK_SIZE until 1 + endColumn / BLOCK_SIZE) {
                val jWidth = blockWidth(jBlock)
                val q0 = jBlock * BLOCK_SIZE
                val qStart = max(startColumn, q0)
                val qEnd = min((jBlock + 1) * BLOCK_SIZE, 1 + endColumn)
                val block = blocks[iBlock * blockColumns + jBlock]
                for (p in pStart until pEnd) {
                    var k = (p - p0) * jWidth + qStart - q0
                    for (q in qStart until qEnd) {
                        visitor.visit(p, q, block!![k])
                        ++k
                    }
                }
            }
        }
        return visitor.end()
    }

    /**
     * Get the height of a block.
     *
     * @param blockRow row index (in block sense) of the block
     * @return height (number of rows) of the block
     */
    private fun blockHeight(blockRow: Int): Int {
        return if (blockRow == blockRows - 1) rows - blockRow * BLOCK_SIZE else BLOCK_SIZE
    }

    /**
     * Get the width of a block.
     *
     * @param blockColumn column index (in block sense) of the block
     * @return width (number of columns) of the block
     */
    private fun blockWidth(blockColumn: Int): Int {
        return if (blockColumn == blockColumns - 1) columns - blockColumn * BLOCK_SIZE else BLOCK_SIZE
    }

    companion object {
        /**
         * Block size.
         */
        const val BLOCK_SIZE = 52

        /**
         * Convert a data array from raw layout to blocks layout.
         *
         *
         * Raw layout is the straightforward layout where element at row i and
         * column j is in array element `rawData[i][j]`. Blocks layout
         * is the layout used in [BlockRealMatrix] instances, where the matrix
         * is split in square blocks (except at right and bottom side where blocks may
         * be rectangular to fit matrix size) and each block is stored in a flattened
         * one-dimensional array.
         *
         *
         *
         * This method creates an array in blocks layout from an input array in raw layout.
         * It can be used to provide the array argument of the [ ][.BlockRealMatrix] constructor.
         *
         *
         * @param rawData Data array in raw layout.
         * @return a new data array containing the same entries but in blocks layout.
         * @throws DimensionMismatchException if `rawData` is not rectangular.
         * @see .createBlocksLayout
         * @see .BlockRealMatrix
         */
        @Throws(DimensionMismatchException::class)
        fun toBlocksLayout(rawData: Array<DoubleArray>): Array<DoubleArray?> {
            val rows = rawData.size
            val columns: Int = rawData[0].size
            val blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE
            val blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE

            // safety checks
            for (i in rawData.indices) {
                val length: Int = rawData[i].size
                if (length != columns) {
                    throw DimensionMismatchException(columns, length)
                }
            }

            // convert array
            val blocks = arrayOfNulls<DoubleArray>(blockRows * blockColumns)
            var blockIndex = 0
            for (iBlock in 0 until blockRows) {
                val pStart = iBlock * BLOCK_SIZE
                val pEnd = min(pStart + BLOCK_SIZE, rows)
                val iHeight = pEnd - pStart
                for (jBlock in 0 until blockColumns) {
                    val qStart = jBlock * BLOCK_SIZE
                    val qEnd = min(qStart + BLOCK_SIZE, columns)
                    val jWidth = qEnd - qStart

                    // allocate new block
                    val block = DoubleArray(iHeight * jWidth)
                    blocks[blockIndex] = block

                    // copy data
                    var index = 0
                    for (p in pStart until pEnd) {
                        ArrayUtil.arraycopy(rawData[p], qStart, block, index, jWidth)
                        index += jWidth
                    }
                    ++blockIndex
                }
            }
            return blocks
        }

        /**
         * Create a data array in blocks layout.
         *
         *
         * This method can be used to create the array argument of the [ ][.BlockRealMatrix] constructor.
         *
         *
         * @param rows    Number of rows in the new matrix.
         * @param columns Number of columns in the new matrix.
         * @return a new data array in blocks layout.
         * @see .toBlocksLayout
         * @see .BlockRealMatrix
         */

        fun createBlocksLayout(rows: Int, columns: Int): Array<DoubleArray?> {
            val blockRows = (rows + BLOCK_SIZE - 1) / BLOCK_SIZE
            val blockColumns = (columns + BLOCK_SIZE - 1) / BLOCK_SIZE
            val blocks = arrayOfNulls<DoubleArray>(blockRows * blockColumns)
            var blockIndex = 0
            for (iBlock in 0 until blockRows) {
                val pStart = iBlock * BLOCK_SIZE
                val pEnd = min(pStart + BLOCK_SIZE, rows)
                val iHeight = pEnd - pStart
                for (jBlock in 0 until blockColumns) {
                    val qStart = jBlock * BLOCK_SIZE
                    val qEnd = min(qStart + BLOCK_SIZE, columns)
                    val jWidth = qEnd - qStart
                    blocks[blockIndex] = DoubleArray(iHeight * jWidth)
                    ++blockIndex
                }
            }
            return blocks
        }
    }
}