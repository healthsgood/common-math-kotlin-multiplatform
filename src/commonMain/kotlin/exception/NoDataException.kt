package exception

import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when the required data is missing.
 *
 * @since 2.2
 */
class NoDataException
/**
 * Construct the exception with a specific context.
 *
 * @param specific Contextual information on what caused the exception.
 */
/**
 * Construct the exception.
 */
constructor(specific: Localizable = LocalizedFormats.NO_DATA) : MathIllegalArgumentException(specific)