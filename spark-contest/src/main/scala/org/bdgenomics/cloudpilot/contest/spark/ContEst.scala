package org.bdgenomics.cloudpilot.contest.spark

import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.{ Variant, AlignmentRecord }

class ContEst(val reads: RDD[AlignmentRecord],
              val variants: RDD[Variant]) {

  def alleleFractions(): Unit = {
    variants.foreach {
      v =>
        v.get
    }
  }
}
