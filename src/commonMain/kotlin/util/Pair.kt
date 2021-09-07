/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package util

/**
 * Generic pair.
 * <br></br>
 * Although the instances of this class are immutable, it is impossible
 * to ensure that the references passed to the constructor will not be
 * modified by the caller.
 *
 * @param <K> Key type.
 * @param <V> Value type.
 * @since 3.0
</V></K> */
class Pair<K, V>(k: K, v: V) {
    /**
     * Get the key.
     *
     * @return the key (first element of the pair).
     */
    /**
     * Get the first element of the pair.
     *
     * @return the first element of the pair.
     * @since 3.1
     */
    /**
     * Key.
     */
    val first: K?
    /**
     * Get the value.
     *
     * @return the value (second element of the pair).
     */
    /**
     * Get the second element of the pair.
     *
     * @return the second element of the pair.
     * @since 3.1
     */
    /**
     * Value.
     */
    val second: V?

    /**
     * Compare the specified object with this entry for equality.
     *
     * @param o Object.
     * @return `true` if the given object is also a map entry and
     * the two entries represent the same mapping.
     */
    override fun equals(o: Any?): Boolean {
        if (this === o) {
            return true
        }
        return if (o !is Pair<*, *>) {
            false
        } else {
            val oP = o
            (if (first == null) oP.first == null else first == oP.first) &&
                    if (second == null) oP.second == null else second == oP.second
        }
    }

    /**
     * Compute a hash code.
     *
     * @return the hash code value.
     */
    override fun hashCode(): Int {
        var result = if (first == null) 0 else first.hashCode()
        val h = if (second == null) 0 else second.hashCode()
        result = 37 * result + h xor (h ushr 16)
        return result
    }

    /**
     * {@inheritDoc}
     */
    override fun toString(): String {
        return "[" + first + ", " + second + "]"
    }

    companion object {
        /**
         * Convenience factory method that calls the
         * [constructor][.Pair].
         *
         * @param <K> the key type
         * @param <V> the value type
         * @param k   First element of the pair.
         * @param v   Second element of the pair.
         * @return a new `Pair` containing `k` and `v`.
         * @since 3.3
        </V></K> */
        fun <K, V> create(k: K, v: V): Pair<K, V> {
            return Pair(k, v)
        }
    }

    /**
     * Create an entry representing a mapping from the specified key to the
     * specified value.
     *
     * @param k Key (first element of the pair).
     * @param v Value (second element of the pair).
     */
    init {
        first = k
        second = v
    }
}