package org.bdgenomics.cloudpilot.contest.spark

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{ VariantContext, SnpTable }
import org.bdgenomics.avocado.models.{ Observation, ReadObservation }
import org.bdgenomics.formats.avro.AlignmentRecord

import org.bdgenomics.avocado.genotyping.MutectGenotyper

import org.apache.spark.SparkContext._
import org.bdgenomics.adam.rdd.ADAMContext._

class SomaticCallingPipeline(val tumorReads: RDD[AlignmentRecord],
                             val normalReads: RDD[AlignmentRecord],
                             val snpTable: SnpTable,
                             val homozygousVariantSites: RDD[VariantSite])(
                                 implicit sc: SparkContext, mutect: MutectGenotyper) {

  val tumorDedup = tumorReads.adamMarkDuplicates()
  val normalDedup = normalReads.adamMarkDuplicates()

  val taggedDedup = tumorDedup.map(tagTumor).union(normalDedup.map(tagNormal))
  val taggedRealigned = taggedDedup.adamRealignIndels()

  val tumorRealigned = taggedRealigned.filter(isTumor)
  val normalRealigned = taggedRealigned.filter(isNormal)

  val contest = new ContEst(tumorRealigned, homozygousVariantSites)

  val broadcastSnpTable = sc.broadcast(snpTable)
  val tumorCorrected = tumorRealigned.adamBQSR(broadcastSnpTable)
  val normalCorrected = normalRealigned.adamBQSR(broadcastSnpTable)

  /**
   * Runs the Mutect algorithm on the tumor and normal reads and inferring the locations of
   * somatic point mutations.
   *
   * @return An RDD of VariantContexts, each of which represents a called somatic mutation.
   */
  def callVariants(): RDD[VariantContext] = {
    val tumorObservations: RDD[Observation] = tumorCorrected.map(rec => new ReadObservation(rec))
    val normalObservations: RDD[Observation] = normalCorrected.map(rec => new ReadObservation(rec))
    val observations: RDD[Observation] = tumorObservations.union(normalObservations)

    mutect.genotype(observations)
  }

  def unitIntervalBins(numBins: Int): Seq[Double] = {
    val binWidth = 1.0 / numBins
    (0 until (numBins + 1)).map(i => i * binWidth)
  }

  /**
   * Estimates, using the ContEst algorithm, the amount of contamination in the
   * tumor reads.
   *
   * @param numBins Determines the number of bins into which the unit interval is split,
   *                defining the points at which the likelihood of the contamination
   *                value 'c' is estimated and therefore the accuracy of that estimate.
   * @return The estimate of c
   */
  def estimateContamination(numBins: Int): Double = {
    val contest = new ContEst(tumorReads, homozygousVariantSites)
    val bins = unitIntervalBins(numBins)

    val lls: Seq[(Double, Double)] =
      bins.map(c => (c, contest.logLikelihood(c))).sortBy(_._2).reverse

    // TODO: calculate the empirical confidence interval
    lls.head._1
  }

  private def tagRead(newTag: String)(read: AlignmentRecord): AlignmentRecord = {
    read.setAttributes(read.getAttributes + " %s".format(newTag))
    read
  }

  private def tagNormal(r: AlignmentRecord): AlignmentRecord = tagRead("XN")(r)
  private def tagTumor(r: AlignmentRecord): AlignmentRecord = tagRead("XT")(r)

  private def isNormal(read: AlignmentRecord): Boolean =
    read.tags.exists(_.tag == "XN")

  private def isTumor(read: AlignmentRecord): Boolean =
    read.tags.exists(_.tag == "XT")

}
