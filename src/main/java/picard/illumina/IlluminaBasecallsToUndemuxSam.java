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

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.illumina.parser.IlluminaFileUtil;
import picard.illumina.parser.IlluminaXMLParser;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.AdapterPair;
import picard.util.IlluminaUtil.IlluminaAdapterPair;

import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * IlluminaBasecallsToSam transforms a lane of Illumina data file formats (bcl, locs, clocs, qseqs, etc.) into
 * SAM or BAM file format.
 * <p/>
 * In this application, barcode data is read from Illumina data file groups, each of which is associated with a tile.
 * Each tile may contain data for any number of barcodes, and a single barcode's data may span multiple tiles.  Once the
 * barcode data is collected from files, each barcode's data is written to its own SAM/BAM.  The barcode data must be
 * written in order; this means that barcode data from each tile is sorted before it is written to file, and that if a
 * barcode's data does span multiple tiles, data collected from each tile must be written in the order of the tiles
 * themselves.
 * <p/>
 * This class employs a number of private subclasses to achieve this goal.  The TileReadAggregator controls the flow
 * of operation.  It is fed a number of Tiles which it uses to spawn TileReaders.  TileReaders are responsible for
 * reading Illumina data for their respective tiles from disk, and as they collect that data, it is fed back into the
 * TileReadAggregator.  When a TileReader completes a tile, it notifies the TileReadAggregator, which reviews what was
 * read and conditionally queues its writing to disk, baring in mind the requirements of write-order described in the
 * previous paragraph.  As writes complete, the TileReadAggregator re-evaluates the state of reads/writes and may queue
 * more writes.  When all barcodes for all tiles have been written, the TileReadAggregator shuts down.
 * <p/>
 * The TileReadAggregator controls task execution using a specialized ThreadPoolExecutor.  It accepts special Runnables
 * of type PriorityRunnable which allow a priority to be assigned to the runnable.  When the ThreadPoolExecutor is
 * assigning threads, it gives priority to those PriorityRunnables with higher priority values.  In this application,
 * TileReaders are assigned lowest priority, and write tasks are assigned high priority.  It is designed in this fashion
 * to minimize the amount of time data must remain in memory (write the data as soon as possible, then discard it from
 * memory) while maximizing CPU usage.
 *
 * @author jburke@broadinstitute.org
 * @author mccowan@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = IlluminaBasecallsToUndemuxSam.USAGE_SUMMARY + IlluminaBasecallsToUndemuxSam.USAGE_DETAILS,
        oneLineSummary = IlluminaBasecallsToUndemuxSam.USAGE_SUMMARY,
        programGroup = BaseCallingProgramGroup.class
)
@DocumentedFeature
public class IlluminaBasecallsToUndemuxSam extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Transforms raw Illumina sequencing data into an unmapped, undemultiplexed SAM or BAM file.";
    static final String USAGE_DETAILS = "<p>The IlluminaBasecallsToUndemuxSam program collects reads across all " +
            "of the tiles of a lane to produce an unmapped SAM/BAM file.  An unmapped BAM file is often referred to as a uBAM.</p> " +

            "<h4>Usage example:</h4>" +
            "<pre>" +
            "" +
            "java -jar picard.jar IlluminaBasecallsToSam \\<br />" +
            "      BASECALLS_DIR=/BaseCalls/ \\<br />" +
            "      LANE=001 \\<br />" +
            "      READ_STRUCTURE=25T8B25T \\<br />" +
            "      RUN_BARCODE=run15 \\<br />" +
            "      IGNORE_UNEXPECTED_BARCODES=true \\<br />" +
            "      LIBRARY_PARAMS=library.params " +
            "</pre>" +
            "<hr />";

    @Argument(doc = "The run folder. ", shortName = "RD")
    public File RUN_DIR;

    @Argument(doc = "Lane number. ", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Argument(doc = "The output SAM or BAM file. Format is determined by extension.",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "The name of the sequencing center that produced the reads.  Used to set the RG.CN tag.", optional = true)
    public String SEQUENCING_CENTER = "BI";

    @Argument(doc = ReadStructure.PARAMETER_DOC, shortName = "RS", optional = true)
    public String READ_STRUCTURE;

    @Argument(doc = "Which adapters to look for in the read.")
    public List<IlluminaAdapterPair> ADAPTERS_TO_CHECK = new ArrayList<>(
            Arrays.asList(IlluminaAdapterPair.INDEXED,
                    IlluminaAdapterPair.DUAL_INDEXED,
                    IlluminaAdapterPair.NEXTERA_V2,
                    IlluminaAdapterPair.FLUIDIGM));

    @Argument(doc = "For specifying adapters other than standard Illumina", optional = true)
    public String FIVE_PRIME_ADAPTER;

    @Argument(doc = "For specifying adapters other than standard Illumina", optional = true)
    public String THREE_PRIME_ADAPTER;

    @Argument(doc = "The number of threads to run in parallel. If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0, then the number of cores used will" +
            " be the number available on the machine less NUM_PROCESSORS.")
    public Integer NUM_PROCESSORS = 0;

    @Argument(doc = "If set, this is the first tile to be processed (used for debugging).  Note that tiles are not processed" +
            " in numerical order.",
            optional = true)
    public Integer FIRST_TILE;

    @Argument(doc = "If set, process no more than this many tiles (used for debugging).", optional = true)
    public Integer TILE_LIMIT;

    @Argument(doc = "If true, call System.gc() periodically.  This is useful in cases in which the -Xmx value passed " +
            "is larger than the available memory.")
    public Boolean FORCE_GC = true;

    @Argument(doc = "Apply EAMSS filtering to identify inappropriately quality scored bases towards the ends of reads" +
            " and convert their quality scores to Q2.")
    public boolean APPLY_EAMSS_FILTER = true;

    @Argument(doc = "Configure SortingCollections to store this many records before spilling to disk. For an indexed" +
            " run, each SortingCollection gets this value/number of indices.")
    public int MAX_READS_IN_RAM_PER_TILE = 1200000;

    @Argument(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown." +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Argument(doc = "Whether to include non-PF reads", shortName = "NONPF", optional = true)
    public boolean INCLUDE_NON_PF_READS = true;

    @Argument(doc = "The tag to use to store any molecular indexes.  If more than one molecular index is found, they will be concatenated and stored here.", optional = true)
    public String MOLECULAR_INDEX_TAG = "RX";

    @Argument(doc = "The tag to use to store any molecular index base qualities.  If more than one molecular index is found, their qualities will be concatenated and stored here " +
            "(.i.e. the number of \"M\" operators in the READ_STRUCTURE)", optional = true)
    public String MOLECULAR_INDEX_BASE_QUALITY_TAG = "QX";

    @Argument(doc = "The list of tags to store each molecular index.  The number of tags should match the number of molecular indexes.", optional = true)
    public List<String> TAG_PER_MOLECULAR_INDEX;

    public String PLATFORM = "ILLUMINA";

    private final Map<String, IlluminaBasecallsToSam.SAMFileWriterWrapper> barcodeSamWriterMap = new HashMap<>();
    private ReadStructure readStructure;
    private BasecallsConverter<IlluminaBasecallsToSam.SAMRecordsForCluster> basecallsConverter;
    private static final Log log = Log.getInstance(IlluminaBasecallsToUndemuxSam.class);
    private File baseCallsDir;

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access args.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     * to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> messages = new ArrayList<>();

        if (!TAG_PER_MOLECULAR_INDEX.isEmpty() && TAG_PER_MOLECULAR_INDEX.size() != readStructure.molecularBarcode.length()) {
            messages.add("The number of tags given in TAG_PER_MOLECULAR_INDEX does not match the number of molecular indexes in READ_STRUCTURE");
        }

        if ((FIVE_PRIME_ADAPTER == null) != (THREE_PRIME_ADAPTER == null)) {
            messages.add("THREE_PRIME_ADAPTER and FIVE_PRIME_ADAPTER must either both be null or both be set.");
        }

        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }


    @Override
    protected int doWork() {
        initialize();
        basecallsConverter.doTileProcessing();
        return 0;
    }

    /**
     * Prepares loggers, initiates garbage collection thread, parses arguments and initialized variables appropriately/
     */
    private void initialize() {
        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);

        IOUtil.assertDirectoryIsReadable(RUN_DIR);
        IOUtil.assertFileIsWritable(OUTPUT);

        baseCallsDir = new File(RUN_DIR,"Data"+File.separator+"Intensities"+File.separator+"BaseCalls");
        IOUtil.assertDirectoryIsReadable(baseCallsDir);

        IlluminaXMLParser illuminaXMLParser = new IlluminaXMLParser(RUN_DIR);

        if (READ_STRUCTURE!=null) {
            readStructure = new ReadStructure(READ_STRUCTURE);
        }else{
            readStructure = illuminaXMLParser.getReadStructure();
        }

        barcodeSamWriterMap.put(null, buildSamFileWriter(OUTPUT, illuminaXMLParser));

        final int numOutputRecords = readStructure.templates.length();

        // Combine any adapters and custom adapter pairs from the command line into an array for use in clipping
        final List<AdapterPair> adapters = new ArrayList<>();
        adapters.addAll(ADAPTERS_TO_CHECK);
        if (FIVE_PRIME_ADAPTER != null && THREE_PRIME_ADAPTER != null) {
            adapters.add(new CustomAdapterPair(FIVE_PRIME_ADAPTER, THREE_PRIME_ADAPTER));
        }

        if (IlluminaFileUtil.hasCbcls(baseCallsDir, LANE)) {
            basecallsConverter = new NewIlluminaBasecallsConverter<>(baseCallsDir, baseCallsDir, LANE, readStructure,
                    barcodeSamWriterMap, false, Math.max(1, MAX_READS_IN_RAM_PER_TILE / numOutputRecords),
                    TMP_DIR, NUM_PROCESSORS,
                    FIRST_TILE, TILE_LIMIT, new IlluminaBasecallsToSam.QueryNameComparator(),
                    new IlluminaBasecallsToSam.Codec(numOutputRecords),
                    IlluminaBasecallsToSam.SAMRecordsForCluster.class, bclQualityEvaluationStrategy, true);
        } else {
            basecallsConverter = new IlluminaBasecallsConverter<>(baseCallsDir, baseCallsDir, LANE, readStructure,
                    barcodeSamWriterMap, false, MAX_READS_IN_RAM_PER_TILE / numOutputRecords, TMP_DIR, NUM_PROCESSORS, FORCE_GC,
                    FIRST_TILE, TILE_LIMIT, new IlluminaBasecallsToSam.QueryNameComparator(), new IlluminaBasecallsToSam.Codec(numOutputRecords),IlluminaBasecallsToSam.SAMRecordsForCluster.class,
                    bclQualityEvaluationStrategy, APPLY_EAMSS_FILTER, INCLUDE_NON_PF_READS, true);
        }
        /*
         * Be sure to pass the outputReadStructure to ClusterDataToSamConverter, which reflects the structure of the output cluster
         * data which may be different from the input read structure (specifically if there are skips).
         */
        final ClusterDataToSamConverter converter = new ClusterDataToUndemuxSamConverter(
                getRunBarcode(illuminaXMLParser),getReadGroupId(illuminaXMLParser),
                basecallsConverter.getFactory().getOutputReadStructure(), adapters)
                .withMolecularIndexTag(MOLECULAR_INDEX_TAG)
                .withMolecularIndexQualityTag(MOLECULAR_INDEX_BASE_QUALITY_TAG)
                .withTagPerMolecularIndex(TAG_PER_MOLECULAR_INDEX);

        basecallsConverter.setConverter(converter);
        log.info("DONE_READING STRUCTURE IS " + readStructure.toString());
    }


    private String getRunBarcode(IlluminaXMLParser illuminaXMLParser){
        return "HWI-"+illuminaXMLParser.getInstrumentId()+"_"+illuminaXMLParser.getExperimentName();
    }


    private String getReadGroupId(IlluminaXMLParser illuminaXMLParser){
        return illuminaXMLParser.getFlowcellId()+"_"+LANE;
    }

    /**
     * Create the list of headers that will be added to the SAMFileHeader for a library with the given sampleBarcodes (or
     * the entire run if sampleBarcodes == NULL).  Note that any value that is null will NOT be added via buildSamFileWriter
     * but is placed in the map in order to be able to query the tags that we automatically add.
     *
     * @param illuminaXMLParser The illumina XML parser used to parse RunParameters.xml and RunInfo.xml files
     * @return A Map of ReadGroupHeaderTags -> Values
     */
    private Map<String, String> buildSamHeaderParameters(IlluminaXMLParser illuminaXMLParser) {
        final Map<String, String> params = new LinkedHashMap<>();

        params.put("PL", PLATFORM);
        params.put("PU", illuminaXMLParser.getIlluminaUniqueRunId());
        params.put("CN", SEQUENCING_CENTER);

        SimpleDateFormat sdf = new SimpleDateFormat("yyMMdd");
        try {
            params.put("DT", new Iso8601Date(sdf.parse(illuminaXMLParser.getIlluminaUniqueRunId().split("_")[0])).toString());
        }catch (ParseException e){
            throw new PicardException("Cannot parse the date in the illumina Run ID "+illuminaXMLParser.getIlluminaUniqueRunId());
        }

        return params;
    }

    /**
     * Create a Program Line for this progrm that will be added to the SAMFileHeader
     *
     * @param previousProgramId The ID of the previous program, generally the Illumina machine RTA software
     * @return A SAMProgramRecord containing this program information
     */
    private SAMProgramRecord buildThisProgramHeaderInfo(String previousProgramId){
        SAMProgramRecord result = new SAMProgramRecord(this.getClass().getSimpleName());
        result.setProgramName(this.getClass().getSimpleName());
        result.setProgramVersion(this.getVersion());
        result.setCommandLine(this.getCommandLine());
        result.setAttribute("DS","Convert Illumina BCL to unaligned BAM or SAM file");

        result.setPreviousProgramGroupId(previousProgramId);
        return result;
    }

    /**
     * Build a SamFileWriter that will write its contents to the output file.
     *
     * @param output           The file to which to write
     * @param illuminaXMLParser The illumina XML parser used to parse RunParameters.xml and RunInfo.xml files
     * @return A SAMFileWriter
     */
    private IlluminaBasecallsToSam.SAMFileWriterWrapper buildSamFileWriter(final File output, IlluminaXMLParser illuminaXMLParser) {
        IOUtil.assertFileIsWritable(output);

        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);

        final SAMProgramRecord thisProgramRecord = buildThisProgramHeaderInfo(illuminaXMLParser.getRtaProgram().getProgramGroupId());

        String rgId = getReadGroupId(illuminaXMLParser);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(rgId);

        rg.setSample(rgId);
        rg.setLibrary(rgId);

        for (final Map.Entry<String, String> tagNameToValue : buildSamHeaderParameters(illuminaXMLParser).entrySet()) {
            if (tagNameToValue.getValue() != null) {
                rg.setAttribute(tagNameToValue.getKey(), tagNameToValue.getValue());
            }
        }

        rg.setProgramGroup(thisProgramRecord.getProgramGroupId());
        header.addReadGroup(rg);

        header.addProgramRecord(illuminaXMLParser.getInstrumentProgram());
        header.addProgramRecord(illuminaXMLParser.getRtaProgram());
        header.addProgramRecord(thisProgramRecord);

        String comment = "READ_STRUCTURE="+readStructure.toString();
        header.addComment(comment);

        return new IlluminaBasecallsToSam.SAMFileWriterWrapper(new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, output));
    }

    public static void main(final String[] args) {
        System.exit(new IlluminaBasecallsToUndemuxSam().instanceMain(args));
    }



}
