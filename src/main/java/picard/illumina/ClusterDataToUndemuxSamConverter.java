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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.ArrayUtils;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.ReadStructure;
import picard.util.AdapterPair;

import java.util.List;

/**
 * Takes ClusterData provided by an IlluminaDataProvider into one or two SAMRecords,
 * as appropriate, and optionally marking adapter sequence.  There is one converter per
 * IlluminaBasecallsToSam run, and all the TileProcessors use the same converter.
 * 
 * @author jburke@broadinstitute.org
 */
public class ClusterDataToUndemuxSamConverter extends ClusterDataToSamConverter {


    /**
     * Constructor
     *
     * @param runBarcode        Used to construct read names.
     * @param readGroupId       If non-null, set RG attribute on SAMRecord to this.
     * @param readStructure     The expected structure (number of reads and indexes,
     *                          and their length) in the read.
     * @param adapters          The list of adapters to check for in the read
     */
    public ClusterDataToUndemuxSamConverter(final String runBarcode,
                                            final String readGroupId,
                                            final ReadStructure readStructure,
                                            final List<AdapterPair> adapters) {
        super(runBarcode,readGroupId,readStructure,adapters, PopulateBarcode.ALWAYS, true);
    }

    /**
     * Creates the SAMRecord for each read in the cluster
     */
    @Override
    public IlluminaBasecallsToSam.SAMRecordsForCluster convertClusterToOutputRecord(final ClusterData cluster) {

        IlluminaBasecallsToSam.SAMRecordsForCluster result = super.convertClusterToOutputRecord(cluster);

        String barcode = null;
        String barcodeQuality = null;

        if (this.hasSampleBarcode && cluster.getMatchedBarcode() == null) {
            byte[] barcodeArray = new byte[0];
            byte[] barcodeQualArray = new byte[0];
            for (int i = 0; i < sampleBarcodeIndices.length; i++) {
                barcodeArray = ArrayUtils.addAll(barcodeArray, cluster.getRead(sampleBarcodeIndices[i]).getBases());
                barcodeQualArray = ArrayUtils.addAll(barcodeQualArray, cluster.getRead(sampleBarcodeIndices[i]).getQualities());

            }
            barcode = StringUtil.bytesToString(barcodeArray).replace('.', 'N');
            barcodeQuality = SAMUtils.phredToFastq(barcodeQualArray);
        }

        if (barcode != null && barcodeQuality != null) {
            boolean first = true;
            for (SAMRecord record : result.records) {
                if (first) {
                    record.setAttribute(SAMTag.BC.name(), barcode);
                    record.setAttribute(SAMTag.QT.name(), barcodeQuality);
                    first = false;
                }else {
                    record.setAttribute(SAMTag.BC.name(),null);
                }
            }
        }

        return result;

    }
}
