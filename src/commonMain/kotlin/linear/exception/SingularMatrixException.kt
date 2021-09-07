package linear.exception

import exception.MathIllegalArgumentException
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when a non-singular matrix is expected.
 *
 * @author haokangkang
 * @since 3.0
 */
class SingularMatrixException : MathIllegalArgumentException(LocalizedFormats.SINGULAR_MATRIX)