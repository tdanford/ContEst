package org.bdgenomics.cloudpilot.contest.spark

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ShuffleRegionJoin
import org.apache.spark.SparkContext._

import org.bdgenomics.formats.avro.{ Variant, AlignmentRecord }

class ContEst(val reads: RDD[AlignmentRecord],
              val variants: RDD[VariantSite]) {

  import ContEst._

  def logLikelihood(c: Double): Double = {
    val joined: RDD[(VariantSite, AlignmentRecord)] =
      ShuffleRegionJoin.partitionAndJoin(variants.keyBy(ContEst.variantRegion), reads.keyBy(ReferenceRegion(_)))

    val readlikelihoods: RDD[(VariantSite, Double)] = joined.map(carry(calculateLogLikelihood(c)))
    val siteLikelihoods: RDD[(VariantSite, Double)] = readlikelihoods.reduceByKey(_ + _)

    siteLikelihoods.map(_._2).reduce(_ + _)
  }

  def calculateLogLikelihood(c: Double)(vs: VariantSite, read: AlignmentRecord): Double =
    vs.logLikelihood(c, read)

}

object ContEst extends Serializable {

  def carry[U, T, V](f: (U, T) => V): ((U, T)) => (U, V) = {
    def carried(p: (U, T)): (U, V) = {
      p match {
        case (u: U, _) => (u, f(p._1, p._2))
      }
    }
    carried
  }

  def lift[U, T, V](f: T => V): ((U, T)) => (U, V) = {
    def lifted(p: (U, T)): (U, V) = {
      p match {
        case (u: U, t: T) =>
          (u, f(t))
      }
    }
    lifted
  }

  def phredToError(phred: Char): Double = {
    val phredValue: Int = phred - '!'
    Math.exp(-10.0 * phredValue)
  }

  def readRegion(read: AlignmentRecord): ReferenceRegion = ReferenceRegion(read)
  def variantRegion(vs: VariantSite): ReferenceRegion =
    ReferenceRegion(vs.variant.getContig.getContigName, vs.variant.getStart, vs.variant.getEnd)
}

case class VariantSite(variant: Variant, populationAlternateFrequency: Double) extends Serializable {

  def logLikelihood(c: Double, read: AlignmentRecord): Double = {
    val offset: Int = (variant.getStart - read.getStart).toInt
    val readBase: String = read.getSequence.substring(offset, offset + 1)
    val error: Double = ContEst.phredToError(read.getQual.charAt(offset))
    val notError = 1.0 - error
    val error3 = error / 3.0
    val f: Double = populationAlternateFrequency

    if (readBase.equals(variant.getReferenceAllele)) {
      val not_c: Double = (1.0 - c) * notError
      val _c: Double = c * (f * error3 + (1.0 - f) * notError)
      not_c + _c

    } else if (readBase.equals(variant.getAlternateAllele)) {
      val not_c: Double = (1.0 - c) * error3
      val _c: Double = c * (f * error3 + (1.0 - f) * notError)
      not_c + _c

    } else {
      error3
    }
  }
}
