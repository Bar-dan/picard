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
import org.apache.commons.lang3.ArrayUtils;
import picard.PicardException;
import picard.illumina.parser.*;
import picard.util.AdapterMarker;
import picard.util.AdapterPair;
import picard.util.IlluminaUtil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

/**
 * Takes ClusterData provided by an IlluminaDataProvider into one or two SAMRecords,
 * as appropriate, and optionally marking adapter sequence.  There is one converter per
 * IlluminaBasecallsToSam run, and all the TileProcessors use the same converter.
 *
 * @author jburke@broadinstitute.org
 */
public class SamToSamMapper {
    private final SamRecordFilter filters = new SolexaNoiseFilter();

    private final boolean hasSampleBarcode;
    private final boolean hasMolecularBarcode;
    private final boolean hasCellBarcode;

    private final boolean isPairedEnd;

    private final ReadStructure readStructure;

    private final AdapterMarker adapterMarker;

    private String MOLECULAR_INDEX_TAG = "RX";
    private String MOLECULAR_INDEX_QUALITY_TAG = "QX";
    private String BARCODE_QUALITY_TAG_NAME = "QT";
    private String BARCODE_TAG_NAME = "BC";
    private String CELL_INDEX_TAG = "CB";
    private String CELL_INDEX_QUALITY_TAG = "CY";

    private static final String MOLECULAR_INDEX_DELIMITER = "-";
    private static final String MOLECULAR_INDEX_QUALITY_DELIMITER = "~";

    private static final String CELL_INDEX_DELIMITER = "-";
    private static final String CELL_INDEX_QUALITY_DELIMITER = "~";

    private static final Character MISSING_BARCODE = '.';
    private static final Character MISSING_BARCODE_BASE = 'N';

    private List<String> tagPerMolecularIndex = Collections.emptyList();

    /**
     * Constructor
     *
     * @param readStructure             The expected structure (number of reads and indexes,
     *                                  and their length) in the read.
     * @param adapters                  The list of adapters to check for in the read
     */
    public SamToSamMapper(final ReadStructure readStructure,
                          final List<AdapterPair> adapters) throws PicardException {

        this.hasSampleBarcode = !readStructure.sampleBarcodes.isEmpty();
        this.hasMolecularBarcode = !readStructure.molecularBarcode.isEmpty();
        this.hasCellBarcode = !readStructure.cellBarcode.isEmpty();

        if (adapters.isEmpty()) {
            this.adapterMarker = null;
        } else {
            this.adapterMarker = new AdapterMarker(adapters.toArray(new AdapterPair[adapters.size()]));
        }

        this.readStructure = readStructure;
        this.isPairedEnd = this.readStructure.templates.length() == 2;

        if (this.readStructure.templates.length()!=1 && this.readStructure.templates.length()!=2) {
            throw new PicardException("There can be only one or two templates in a read structure");
        }
    }

    /**
     * Sets the SAM tag to use to store the molecular index bases.  If multiple molecular indexes exist, it will concatenate them
     * and store them in this tag.
     */
    public SamToSamMapper withMolecularIndexTag(final String molecularIndexTag) {
        if (molecularIndexTag == null) throw new IllegalArgumentException("Molecular index tag was null");
        this.MOLECULAR_INDEX_TAG = molecularIndexTag;
        return this;
    }

    /**
     * Sets the SAM tag to use to store the molecular index base qualities.  If multiple molecular indexes exist, it will concatenate them
     * and store them in this tag.
     */
    public SamToSamMapper withMolecularIndexQualityTag(final String molecularIndexQualityTag) {
        if (molecularIndexQualityTag == null) {
            throw new IllegalArgumentException("Molecular index quality tag was null");
        }
        this.MOLECULAR_INDEX_QUALITY_TAG = molecularIndexQualityTag;
        return this;
    }

    /**
     * Sets the SAM tags to use to store the bases each molecular index.  This will only be used if there are more than one molecular
     * index. If fewer tags are given than molecular indexes found, then the remaining molecular indexes will be concatenated and stored
     * in the last tag.  If more tags are provided than molecular indexes found, the additional tags will not be used.
     */
    public SamToSamMapper withTagPerMolecularIndex(final List<String> tagPerMolecularIndex) {
        if (tagPerMolecularIndex == null) {
            throw new IllegalArgumentException("Null given for tagPerMolecularIndex");
        }
        this.tagPerMolecularIndex = tagPerMolecularIndex;
        return this;
    }

    public SamToSamMapper withBarcodeTag(final String barcodeTag) {
        if (barcodeTag == null) {
            throw new IllegalArgumentException("Null given for barcodeTag");
        }
        this.BARCODE_TAG_NAME = barcodeTag;
        return this;
    }

    public SamToSamMapper withBarcodeQualityTag(final String barcodeQualityTag) {
        if (barcodeQualityTag == null) {
            throw new IllegalArgumentException("Null given for barcodeQualityTag");
        }
        this.BARCODE_QUALITY_TAG_NAME = barcodeQualityTag;
        return this;
    }

    public SamToSamMapper withCellIndexTag(final String cellIndexTag) {
        if (cellIndexTag == null) {
            throw new IllegalArgumentException("Null given for cellIndexTag");
        }
        this.CELL_INDEX_TAG = cellIndexTag;
        return this;
    }

    public SamToSamMapper withCellIndexQualityTag(final String cellIndexQualityTag) {
        if (cellIndexQualityTag == null) {
            throw new IllegalArgumentException("Null given for cellIndexQualityTag");
        }
        this.CELL_INDEX_QUALITY_TAG = cellIndexQualityTag;
        return this;
    }

    /**
     * Creates a new SAM record from the basecall data
     */
    private SAMRecord createSamRecord(final byte[] bases, final byte[] qualities,
                                      final String readName,
                                      final String readGroupId,
                                      final boolean failsVendorQc,
                                      final boolean firstOfPair,
                                      final String unmatchedBarcode,
                                      final String barcodeQuality,
                                      final List<String> molecularIndexes,
                                      final List<String> molecularIndexQualities,
                                      final List<String> cellIndexes,
                                      final List<String> cellIndexQualities) {
        final SAMRecord sam = new SAMRecord(null);
        sam.setReadName(readName);
        sam.setReadBases(bases);
        sam.setBaseQualities(qualities);

        // Flag values
        sam.setReadPairedFlag(isPairedEnd);
        sam.setReadUnmappedFlag(true);
        sam.setReadFailsVendorQualityCheckFlag(failsVendorQc);
        if (isPairedEnd) {
            sam.setMateUnmappedFlag(true);
            sam.setFirstOfPairFlag(firstOfPair);
            sam.setSecondOfPairFlag(!firstOfPair);
        }

        if (filters.filterOut(sam)) {
            sam.setAttribute(ReservedTagConstants.XN, 1);
        }

        sam.setAttribute(SAMTag.RG.name(), readGroupId);

        if (unmatchedBarcode != null) {
            sam.setAttribute(SAMTag.BC.name(), unmatchedBarcode);
            if (barcodeQuality != null ) {
                sam.setAttribute(SAMTag.QT.name(), barcodeQuality);
            }
        }

        if (!molecularIndexes.isEmpty()) {
            if (!this.MOLECULAR_INDEX_TAG.isEmpty()) {
                sam.setAttribute(this.MOLECULAR_INDEX_TAG, String.join(MOLECULAR_INDEX_DELIMITER, molecularIndexes));
            }
            if (!this.MOLECULAR_INDEX_QUALITY_TAG.isEmpty()) {
                sam.setAttribute(this.MOLECULAR_INDEX_QUALITY_TAG, String.join(MOLECULAR_INDEX_QUALITY_DELIMITER, molecularIndexQualities));
            }
            if (!this.tagPerMolecularIndex.isEmpty()) {
                if (tagPerMolecularIndex.size() != molecularIndexes.size()) {
                    throw new PicardException("Found " + molecularIndexes.size() + " molecular indexes but only " + tagPerMolecularIndex.size() + " SAM tags given.");
                }
                for (int i = 0; i < this.tagPerMolecularIndex.size(); i++) {
                    sam.setAttribute(this.tagPerMolecularIndex.get(i), molecularIndexes.get(i));
                }
            }
        }

        if (!cellIndexes.isEmpty()) {
            if (!this.CELL_INDEX_TAG.isEmpty()) {
                sam.setAttribute(this.CELL_INDEX_TAG, String.join(CELL_INDEX_DELIMITER, cellIndexes));
            }
            if (!this.CELL_INDEX_QUALITY_TAG.isEmpty()) {
                sam.setAttribute(this.CELL_INDEX_QUALITY_TAG, String.join(CELL_INDEX_QUALITY_DELIMITER, cellIndexQualities));
            }
        }

        return sam;
    }

    private int appendRead(byte[] read, byte[] readQual, byte[] dest, byte[] destQual, int cursor) throws PicardException {
        int finalLenght = cursor + read.length;

        if (finalLenght > dest.length) throw new PicardException("Only "+dest.length+" bases specified in " +
                "template, but found at least" + finalLenght );
        System.arraycopy(read,0,dest,cursor,read.length);
        System.arraycopy(readQual,0,destQual,cursor,read.length);

        return finalLenght;
    }

    /**
     * Creates the SAMRecord for each read in the cluster
     */
    public IlluminaBasecallsToSam.SAMRecordsForCluster convertClusterToOutputRecord(final SAMRecord firstRecord, final  SAMRecord secondRecord) throws PicardException {

        final IlluminaBasecallsToSam.SAMRecordsForCluster ret = new IlluminaBasecallsToSam.SAMRecordsForCluster(readStructure.templates.length());

        byte[] bases = new byte[readStructure.totalCycles];
        byte[] qualities = new byte[readStructure.totalCycles];

        if (firstRecord == null) throw new PicardException("Got null First SAM record");

        int readBases = appendRead(firstRecord.getReadBases(),firstRecord.getBaseQualities(),bases,qualities,0);
        if (firstRecord.hasAttribute(BARCODE_TAG_NAME)){
            if (!firstRecord.hasAttribute(BARCODE_QUALITY_TAG_NAME))
                throw new PicardException("Could not find barcode quality attribute " + BARCODE_QUALITY_TAG_NAME + " in record " + firstRecord.getReadName());

            byte[] bcode = htsjdk.samtools.util.StringUtil.stringToBytes(firstRecord.getStringAttribute(BARCODE_TAG_NAME));
            byte[] bcodeQual = SAMUtils.fastqToPhred(firstRecord.getStringAttribute(BARCODE_QUALITY_TAG_NAME));

            readBases = appendRead(bcode, bcodeQual, bases, qualities, readBases);
        }

        if (secondRecord!=null){
            readBases = appendRead(secondRecord.getReadBases(), secondRecord.getBaseQualities(), bases, qualities, readBases);
        }

        if (readBases != readStructure.totalCycles) throw new PicardException("Expected "+readStructure.totalCycles +" Cycles, but only "+readBases+" bases found in SAM records");

        // Get and transform the unmatched barcode, if any, to store with the reads
        String unmatchedBarcode = null;
        if (hasSampleBarcode) {
            unmatchedBarcode = getUnmatchedBarcode(bases);
        }

        String barcodeQuality = null;
        if (unmatchedBarcode != null) {
            barcodeQuality = getBarcodeQuality(qualities);
        }

        final List<String> molecularIndexes = new ArrayList<>();
        final List<String> molecularIndexQualities = new ArrayList<>();
        if (hasMolecularBarcode) {
            byte[][] molecularBarcodesBases = getSeqs(bases, readStructure.molecularBarcode.getCycleIndexRanges());
            byte[][] molecularBarcodeQualities = getSeqs(qualities, readStructure.molecularBarcode.getCycleIndexRanges());

            for (int i=0; i<readStructure.molecularBarcode.length(); i++){
                molecularIndexes.add(convertMissingToNoCall(new String(molecularBarcodesBases[i])));
                molecularIndexQualities.add(SAMUtils.phredToFastq(molecularBarcodeQualities[i]));
            }
        }

        final List<String> cellIndexes = new ArrayList<>();
        final List<String> cellIndexQualities = new ArrayList<>();
        if (hasCellBarcode) {
            byte[][] cellBarcodeBases = getSeqs(bases, readStructure.cellBarcode.getCycleIndexRanges());
            byte[][] cellBarcodeQualities = getSeqs(qualities, readStructure.cellBarcode.getCycleIndexRanges());

            for (int i=0; i<readStructure.cellBarcode.length(); i++){
                cellIndexes.add(convertMissingToNoCall(new String(cellBarcodeBases[i])));
                cellIndexQualities.add(SAMUtils.phredToFastq(cellBarcodeQualities[i]));
            }
        }

        byte[][] templates = getSeqs(bases, readStructure.templates.getCycleIndexRanges());
        byte[][] templateQualities = getSeqs(qualities, readStructure.templates.getCycleIndexRanges());

        final SAMRecord firstOfPair = createSamRecord(
                templates[0],
                templateQualities[0],
                firstRecord.getReadName(),
                firstRecord.getReadGroup().getReadGroupId(),
                firstRecord.getReadFailsVendorQualityCheckFlag(),
                true,
                unmatchedBarcode,
                barcodeQuality,
                molecularIndexes,
                molecularIndexQualities,
                cellIndexes,
                cellIndexQualities);
        ret.records[0] = firstOfPair;


        SAMRecord secondOfPair = null;

        if (isPairedEnd) {
            secondOfPair = createSamRecord(
                    templates[1],
                    templateQualities[1],
                    firstRecord.getReadName(),
                    firstRecord.getReadGroup().getReadGroupId(),
                    firstRecord.getReadFailsVendorQualityCheckFlag(),
                    true,
                    null,
                    null,
                    new ArrayList<>(),
                    new ArrayList<>(),
                    new ArrayList<>(),
                    new ArrayList<>());
            ret.records[1] = secondOfPair;
        }

        if (adapterMarker != null) {
            // Clip the read
            if (isPairedEnd) {
                adapterMarker.adapterTrimIlluminaPairedReads(firstOfPair, secondOfPair);
            } else {
                adapterMarker.adapterTrimIlluminaSingleRead(firstOfPair);
            }
        }
        return ret;
    }

    private String getBarcodeQuality(byte[] qualities) {

        final StringJoiner barcodeQ = new StringJoiner(MOLECULAR_INDEX_QUALITY_DELIMITER);
        byte[][] barcodeQualities = getSeqs(qualities, readStructure.sampleBarcodes.getCycleIndexRanges());

        for (byte[] barcodeQuality : barcodeQualities) {
            barcodeQ.add(SAMUtils.phredToFastq(barcodeQuality));
        }
        return barcodeQ.toString();
    }

    private String getUnmatchedBarcode(byte[] bases) {
        return convertMissingToNoCall(IlluminaUtil.barcodeSeqsToString(getSeqs(bases, readStructure.sampleBarcodes.getCycleIndexRanges())));
    }

    private byte[][] getSeqs (byte[] array, Range[] cycleRanges){

        final byte[][] result = new byte[cycleRanges.length][];
        for (int i = 0; i < cycleRanges.length; i++) {
            Range range = cycleRanges[i];
            result[i] = ArrayUtils.subarray(array, range.start, range.end+1);
        }
        return result;
    }

    private static String convertMissingToNoCall(final String barcode){
        return barcode.replace(MISSING_BARCODE, MISSING_BARCODE_BASE);
    }
}
