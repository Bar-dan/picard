/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.illumina;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SolexaNoiseFilter;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.PicardException;
import picard.fastq.IlluminaReadNameEncoder;
import picard.fastq.ReadNameEncoder;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.ReadData;
import picard.illumina.parser.ReadStructure;
import picard.util.AdapterMarker;
import picard.util.AdapterPair;
import picard.util.IlluminaUtil;

import java.util.*;

/**
 * Takes ClusterData provided by an IlluminaDataProvider into one or two SAMRecords,
 * as appropriate, and optionally marking adapter sequence.  There is one converter per
 * IlluminaBasecallsToSam run, and all the TileProcessors use the same converter.
 *
 * @author jburke@broadinstitute.org
 */
public class ClusterDataToMultiplexSamConverter implements
        IlluminaBasecallsConverter.ClusterDataConverter<IlluminaBasecallsToSam.SAMRecordsForCluster> {

    private final String readGroupId;
    private final SamRecordFilter filters = new SolexaNoiseFilter();
    private final boolean isPairedEnd;
    private final boolean hasSampleBarcode;
    private final int[] templateIndices;
    private final int[] sampleBarcodeIndices;

    private final int outputRecordsPerCluster;
    private final ReadNameEncoder readNameEncoder;

    private static final Character MISSING_BARCODE = '.';
    private static final Character MISSING_BARCODE_BASE = 'N';

    private static final String BARCODE_QUALITY_DELIMITER = "";
    private static final String BARCODE_DELIMITER = "";

    /**
     * Constructor
     *
     * @param runBarcode                Used to construct read names.
     * @param readGroupId               If non-null, set RG attribute on SAMRecord to this.
     * @param readStructure             The expected structure (number of reads and indexes,
     *                                  and their length) in the read.
     */
    public ClusterDataToMultiplexSamConverter(final String runBarcode,
                                              final String readGroupId,
                                              final ReadStructure readStructure) {

        this.readGroupId = readGroupId;

        this.readNameEncoder = new IlluminaReadNameEncoder(runBarcode);

        this.isPairedEnd = readStructure.templates.length() == 2;
        this.hasSampleBarcode = !readStructure.sampleBarcodes.isEmpty();

        this.templateIndices = readStructure.templates.getIndices();
        this.sampleBarcodeIndices = readStructure.sampleBarcodes.getIndices();
        this.outputRecordsPerCluster = readStructure.templates.length();
    }

    /**
     * Creates a new SAM record from the basecall data
     */
    private SAMRecord createSamRecord(final ReadData readData, final String readName, final boolean isPf, final boolean firstOfPair,
                                      final String unmatchedBarcode, final String barcodeQuality) {

        final SAMRecord sam = new SAMRecord(null);
        sam.setReadName(readName);
        sam.setReadBases(readData.getBases());
        sam.setBaseQualities(readData.getQualities());

        // Flag values
        sam.setReadPairedFlag(isPairedEnd);
        sam.setReadUnmappedFlag(true);
        sam.setReadFailsVendorQualityCheckFlag(!isPf);
        if (isPairedEnd) {
            sam.setMateUnmappedFlag(true);
            sam.setFirstOfPairFlag(firstOfPair);
            sam.setSecondOfPairFlag(!firstOfPair);
        }

        if (filters.filterOut(sam)) {
            sam.setAttribute(ReservedTagConstants.XN, 1);
        }

        if (this.readGroupId != null) {
            sam.setAttribute(SAMTag.RG.name(), readGroupId);
        }

        if (unmatchedBarcode != null) {
            sam.setAttribute(SAMTag.BC.name(), unmatchedBarcode);
            if (barcodeQuality != null ) {
                sam.setAttribute(SAMTag.QT.name(), barcodeQuality);
            }
        }

        return sam;
    }

    /**
     * Creates the SAMRecord for each read in the cluster
     */
    public IlluminaBasecallsToSam.SAMRecordsForCluster convertClusterToOutputRecord(final ClusterData cluster) {

        final IlluminaBasecallsToSam.SAMRecordsForCluster ret = new IlluminaBasecallsToSam.SAMRecordsForCluster(outputRecordsPerCluster);
        final String readName = readNameEncoder.generateReadName(cluster, null); // Use null here to prevent /1 or /2 suffixes on read name.

        // Get and transform the unmatched barcode, if any, to store with the reads
        String unmatchedBarcode = null;
        if (hasSampleBarcode) {
            unmatchedBarcode = getUnmatchedBarcode(cluster);
        }

        String barcodeQuality = null;
        if (unmatchedBarcode != null) {
            barcodeQuality = getBarcodeQuality(cluster);
        }

        final SAMRecord firstOfPair = createSamRecord(
                cluster.getRead(templateIndices[0]),
                readName,
                cluster.isPf(),
                true,
                unmatchedBarcode,
                barcodeQuality);

        ret.records[0] = firstOfPair;

        SAMRecord secondOfPair = null;

        if (isPairedEnd) {
            secondOfPair = createSamRecord(
                    cluster.getRead(templateIndices[1]),
                    readName,
                    cluster.isPf(),
                    false,
                    null,
                    null);
            ret.records[1] = secondOfPair;
        }

        return ret;
    }

    private String getBarcodeQuality(ClusterData cluster) {
        final StringJoiner barcodeQ = new StringJoiner(BARCODE_QUALITY_DELIMITER);

        for (int sampleBarcodeIndex : sampleBarcodeIndices) {
            barcodeQ.add(SAMUtils.phredToFastq(cluster.getRead(sampleBarcodeIndex).getQualities()));
        }
        return barcodeQ.toString();
    }

    private String getUnmatchedBarcode(ClusterData cluster) {
        //return convertMissingToNoCall(IlluminaUtil.barcodeSeqsToString(getBarcodeSeqs(cluster))); //<- This would insert a separator in the bcode
        return convertMissingToNoCall(IlluminaUtil.byteArrayToString(getBarcodeSeqs(cluster), BARCODE_DELIMITER));
    }

    private byte[][] getBarcodeSeqs(ClusterData cluster) {
        final byte[][] barcode = new byte[sampleBarcodeIndices.length][];
        for (int i = 0; i < sampleBarcodeIndices.length; i++) {
            barcode[i] = cluster.getRead(sampleBarcodeIndices[i]).getBases();
        }
        return barcode;
    }

    private static String convertMissingToNoCall(final String barcode){
        return barcode.replace(MISSING_BARCODE, MISSING_BARCODE_BASE);
    }
}