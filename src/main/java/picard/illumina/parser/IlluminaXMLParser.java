package picard.illumina.parser;

import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.Log;
import org.w3c.dom.*;
import org.xml.sax.SAXException;
import picard.PicardException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

public class IlluminaXMLParser {
    private static final Log log = Log.getInstance(IlluminaXMLParser.class);
    private final static String RUN_PARAMETERS_FILE = "RunParameters.xml";
    private final static String RUN_INFO_FILE = "RunInfo.xml";


    private SAMProgramRecord rtaProgram = null;
    private SAMProgramRecord instrumentProgram = null;

    private String flowcellId = null;
    private String instrumentId = null;
    private String illuminaUniqueRunId = null;
    private String experimentName = null;
    private ReadStructure readStructure = null;

    public IlluminaXMLParser(File baseFolder) {
        File runParametersFile = null;
        File runInfoFile = null;

        for (File file : baseFolder.listFiles()) {
            if (file.getName().equalsIgnoreCase(RUN_PARAMETERS_FILE)) {
                runParametersFile = file;
            } else if (file.getName().equalsIgnoreCase(RUN_INFO_FILE)) {
                runInfoFile = file;
            }
        }

        if (runParametersFile == null)
            throw new PicardException("Could not find RunParameters.xml in run folder " + baseFolder.getAbsolutePath());
        if (runInfoFile == null)
            throw new PicardException("Could not find RunInfo.xml in run folder " + baseFolder.getAbsolutePath());


        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        DocumentBuilder db;
        try {
            db = dbf.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            throw new PicardException("Error while initializing XML document builder", e);
        }

        try {
            Document doc = db.parse(runParametersFile);
            parseParameterNode(doc);
        } catch (IOException e) {
            throw new PicardException("Error while reading RunParameters.xml file " + runParametersFile.getAbsolutePath(), e);
        } catch (SAXException e) {
            throw new PicardException("Error while parsing RunParameters.xml file " + runParametersFile.getAbsolutePath(), e);
        }

        try {
            Document doc = db.parse(runInfoFile);
            parseRunInfo(doc);
        } catch (IOException e) {
            throw new PicardException("Error while reading RunInfo.xml file " + runInfoFile.getAbsolutePath(), e);
        } catch (SAXException e) {
            throw new PicardException("Error while parsing RunInfo.xml file " + runInfoFile.getAbsolutePath(), e);
        }

    }

    private Node getNodeByTag(Document document, String... tags) {
        NodeList nodeList;

        for (String tag : tags) {
            nodeList = document.getElementsByTagName(tag);
            if (nodeList != null && (nodeList.getLength() != 0)) {
                return nodeList.item(0);
            }
        }

        return null;
    }

    private void parseParameterNode(Document runParameters) {
        Node applicationNameNode = getNodeByTag(runParameters, "ApplicationName", "Application");
        Node applicationVersionNode = getNodeByTag(runParameters, "ApplicationVersion");
        Node rtaVersionNode = getNodeByTag(runParameters, "RTAVersion", "RtaVersion", "rtaVersion");
        Node experimentNameNode = getNodeByTag(runParameters, "ExperimentName");

        if (applicationNameNode == null)
            throw new PicardException("Could not find application name in RunParameters.xml");
        if (applicationVersionNode == null)
            throw new PicardException("Could not find application version in RunParameters.xml");
        if (rtaVersionNode == null) throw new PicardException("Could not find RTA version in RunParameters.xml");
        if (rtaVersionNode == null) throw new PicardException("Could not find Experiment Name in RunParameters.xml");

        experimentName = experimentNameNode.getFirstChild().getNodeValue();

        String applicationName = applicationNameNode.getFirstChild().getNodeValue();
        String applicationVersion = applicationVersionNode.getFirstChild().getNodeValue();

        instrumentProgram = new SAMProgramRecord("SCS");
        instrumentProgram.setProgramVersion(applicationVersion);
        instrumentProgram.setProgramName(applicationName);
        instrumentProgram.setAttribute("DS", "Controlling software on instrument");

        String rtaVersion = rtaVersionNode.getFirstChild().getNodeValue();

        rtaProgram = new SAMProgramRecord("basecalling");
        rtaProgram.setProgramVersion(rtaVersion);
        rtaProgram.setProgramName("RTA");
        rtaProgram.setAttribute("DS", "Basecalling Package");
        rtaProgram.setPreviousProgramGroupId(instrumentProgram.getId());


    }

    private void parseRunInfo(Document runInfo) {
        final NodeList reads = runInfo.getElementsByTagName("Read");
        final List<ReadDescriptor> descriptors = new ArrayList<>(reads.getLength());
        for (int i = 0; i < reads.getLength(); i++) {
            final Node read = reads.item(i);
            final NamedNodeMap attributes = read.getAttributes();
            final int readNumber = Integer.parseInt(attributes.getNamedItem("Number").getNodeValue());
            final int numCycles = Integer.parseInt(attributes.getNamedItem("NumCycles").getNodeValue());
            final boolean isIndexedRead = attributes.getNamedItem("IsIndexedRead").getNodeValue().toUpperCase().equals("Y");
            if (readNumber != i + 1)
                throw new PicardException("Read number in RunInfo.xml was out of order: " + (i + 1) + " != " + readNumber);
            descriptors.add(new ReadDescriptor(numCycles, isIndexedRead ? ReadType.Barcode : ReadType.Template));
        }
        readStructure = new ReadStructure(descriptors);

        Node illuminaRunNode = getNodeByTag(runInfo, "Run");
        if (illuminaRunNode == null) {
            throw new PicardException("Error while retrieving the run Id from runInfo.xml file");
        }
        illuminaUniqueRunId = illuminaRunNode.getAttributes().getNamedItem("Id").getNodeValue();
        if (illuminaUniqueRunId == null) {
            throw new PicardException("Error while retrieving the run Id from runInfo.xml file");
        }

        Node flowcellIdNode = getNodeByTag(runInfo, "Flowcell");
        if (flowcellIdNode == null) {
            throw new PicardException("Error while retrieving the flowcell Id from runInfo.xml file");
        }
        flowcellId = flowcellIdNode.getFirstChild().getNodeValue();

        Node instrumentIdNode = getNodeByTag(runInfo, "Instrument");
        if (instrumentIdNode == null) {
            throw new PicardException("Error while retrieving the Instrument Id from runInfo.xml file");
        }
        instrumentId = instrumentIdNode.getFirstChild().getNodeValue();


    }


    public SAMProgramRecord getInstrumentProgram() {
        return instrumentProgram;
    }


    public SAMProgramRecord getRtaProgram() {
        return rtaProgram;
    }

    public String getFlowcellId() {
        return flowcellId;
    }

    public String getInstrumentId() {
        return instrumentId;
    }

    public String getIlluminaUniqueRunId() {
        return illuminaUniqueRunId;
    }

    public ReadStructure getReadStructure() {
        return readStructure;
    }

    public String getExperimentName() {
        return experimentName;
    }
}
