package leastsquares

import linear.RealVector

/**
 * Interface for validating a set of model parameters.
 *
 * @author haokangkang
 * @since 3.4
 */
interface ParameterValidator {
    /**
     * Validates the set of parameters.
     *
     * @param params Input parameters.
     * @return the validated values.
     */
    fun validate(params: RealVector?): RealVector?
}