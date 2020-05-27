package BLAST_QC_JAVA;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import org.apache.commons.cli.*;

public class BLASTQC {

    private CommandLineOptions clOptions = new CommandLineOptions();
    private BufferedWriter  hits;
    private BufferedWriter  nohits;
    private BufferedWriter  header;

    public static void main(String[] args) {
        BLASTQC sys = new BLASTQC(args);
    }

    public BLASTQC(String[] args) {
        CLI(args);
        String out = new String();

        try {
            
            hits = new BufferedWriter(new FileWriter(clOptions.getOutfile() + ".hits.txt"));
            nohits = new BufferedWriter(new FileWriter(clOptions.getOutfile() + ".nohits.txt"));
            header = new BufferedWriter(new FileWriter(clOptions.getOutfile() + ".hits.header"));

            if(clOptions.getFileformat().equals("XML")) {

                out = "query_name\tquery_length\taccession_number\tsubject_length\tsubject_description\tE_value"
                + "\tbit_score\tframe\tquery_start\tquery_end\thit_start\thit_end\t%_conserved\t%_identity\n";
                hits.append(out);
                out = "query_name\n";
                nohits.write(out);
                out = "query_name\tsubject_description\n";
                header.write(out);
                parseXML();
            }
            else if(clOptions.getFileformat().equals("tab")) {
                out = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n";
                hits.write(out);
                out = "qseqid\n";
                nohits.write(out);
                out = "qseqid\tsseqid\n";
                header.write(out);
                parseTab();
            }

            hits.close();
            nohits.close();
            header.close();

        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void CLI(String[] args) {
        Options options = new Options();

        Option filename  = Option.builder("f")
                                .hasArg()
                                .required(true)
                                .desc("Specifiy the Blast XML results input file.\n(required)")
                                .build();

        Option fileformat  = Option.builder("ff")
                                .longOpt("fileformat")
                                .hasArg()
                                .required(true)
                                .desc("Specifiy the Blast Results file format. XML or tab.\n(required)")
                                .build();
        
        Option outfile  = Option.builder("o")
                                .longOpt("output")
                                .hasArg()
                                .desc("Specify the output file base name (no extension). Defaults to base name of input file.")
                                .build();

        Option type  = Option.builder("t")
                                .longOpt("type")
                                .hasArg()
                                .required(true)
                                .desc("Specify what type of BLAST you are running\n(Protein or Nucleotide). (required)")
                                .build();

        Option number  = Option.builder("n")
                                .longOpt("number")
                                .hasArg()
                                .optionalArg(true)
                                .desc("Specify the number of hits to return per query sequence. Defaults to return all hits that fit input threshold(s).\n(Int value)")
                                .build();

        Option evalue  = Option.builder("e")
                                .longOpt("evalue")
                                .hasArg()
                                .optionalArg(true)
                                .desc("Specify an evalue threshold.\n(Maximum acceptable evalue) (Float value)")
                                .build();           
                                
        Option definition  = Option.builder("d")
                                .longOpt("definition")
                                .hasArg()
                                .optionalArg(true)
                                .desc("Specify a threshold in the level of definition provided. This is defined by how many separate lines are present in the Hit definition '<Hit_def>' of the XML file.\n(Int value)")
                                .build();           
        
        Option order  = Option.builder("or")
                                .longOpt("order")
                                .hasArg()
                                .optionalArg(true)
                                .desc("Specify the order of the results. By lowest evalue, highest bitscore, highest percent identity or most detailed definition data.\n(default: by evalue 'e')")
                                .build();                                   

        Option erange  = Option.builder("er")
                                .longOpt("erange")
                                .hasArg()
                                .optionalArg(true)
                                .desc("Sets a range of acceptable deviation from the lowest evalue hit in which a more detailed definition would be prefered. Must be ordered by evalue.")
                                .build();   

        Option brange  = Option.builder("br")
                                .longOpt("brange")
                                .hasArg()
                                .optionalArg(true)
                                .desc("Sets a range of acceptable deviation from the highest bitscore hit in which a more detailed definition would be prefered. Must be ordered by bitscore.")
                                .build();  
        
        Option irange  = Option.builder("ir")
                                .longOpt("irange")
                                .hasArg()
                                .optionalArg(true)
                                .desc("Sets a range of acceptable deviation from the highest percent identity hit in which a more detailed definition would be prefered. Must be ordered by percent identity.")
                                .build();  

        options.addOption(filename).addOption(fileformat).addOption(outfile).addOption(type).addOption(number)
        .addOption(evalue).addOption(definition).addOption(order).addOption(erange).addOption(brange).addOption(irange);

            
        if(args.length == 1) {
            System.exit(-1);
        }
       
        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine clo = parser.parse(options, args);

            boolean hasReqArgs = true;
            if(!clo.hasOption("f")) {
                System.out.println("Must specify a BLAST results file using the -f commandline argument.\n");
                hasReqArgs = false;
            }

            if(!clo.hasOption("ff")) {
                System.out.println("Must specify a result file format (XML or tab) using the -ff commandline argument.\n");
                hasReqArgs = false;
            }

            if(!clo.hasOption("t")) {
                System.out.println("Must specify the type of BLAST you are using (n or p) using the -t commandline argument.\n");
                hasReqArgs = false;
            }

            if(!hasReqArgs) {
                System.out.println("Use -h for help menu. Please retry.\n");
                System.exit(2);
            }

            if(clo.hasOption("er") && !(clo.getOptionValue("or", "e").trim().equals("e"))) {
                System.out.println("erange cannot be used. Must order by evalue if this functionality is desired."
                            + "\nuse \"-h\" or \"--help\" to display help menu.");
                System.exit(3);
            }
            if(clo.hasOption("br") && !(clo.getOptionValue("or", "e").trim().equals("b"))) {
                System.out.println("brange cannot be used. Must order by bitscore if this functionality is desired."
                            + "\nuse \"-h\" or \"--help\" to display help menu.");
                System.exit(4);
            }
            if(clo.hasOption("ir") && !(clo.getOptionValue("or", "e").trim().equals("i"))) {
                System.out.println("irange cannot be used. Must order by identity if this functionality is desired."
                            + "\nuse \"-h\" or \"--help\" to display help menu.");
                System.exit(5);
            }

            clOptions.setFilename(clo.getOptionValue("f").trim());
            clOptions.setFileformat(clo.getOptionValue("ff").trim());
            clOptions.setType(clo.getOptionValue("t").trim().charAt(0));

            if(!clo.hasOption("o"))
                clOptions.setOutfile(clOptions.getFilename());
            else 
                clOptions.setOutfile(clo.getOptionValue("o").trim());
            if(!clo.hasOption("n"))
                clOptions.setNumber(0);
            else 
                clOptions.setNumber(Integer.parseInt(clo.getOptionValue("n").trim()));
            if(!clo.hasOption("e"))
                clOptions.setEvalue(Double.MAX_VALUE);
            else
                clOptions.setEvalue(Double.parseDouble(clo.getOptionValue("e").trim()));
            if(!clo.hasOption("b"))
                clOptions.setBitscore(-1);
            else
                clOptions.setBitscore(Double.parseDouble(clo.getOptionValue("b").trim()));
            if(!clo.hasOption("i"))
                clOptions.setIdentity(-1);
            else
                clOptions.setIdentity(Double.parseDouble(clo.getOptionValue("i").trim()));
            if(!clo.hasOption("d"))
                clOptions.setDefinition(-1);
            else
                clOptions.setNumber(Integer.parseInt(clo.getOptionValue("d").trim()));
            if(!clo.hasOption("or"))
                clOptions.setOrder('e');
            else
                clOptions.setOrder(clo.getOptionValue("or").trim().charAt(0));
            if(!clo.hasOption("er"))
                clOptions.setErange(0);
            else
                clOptions.setErange(Double.parseDouble(clo.getOptionValue("er").trim()));
            if(!clo.hasOption("br"))
                clOptions.setBrange(0);
            else
                clOptions.setBrange(Double.parseDouble(clo.getOptionValue("br").trim()));
            if(!clo.hasOption("ir"))
                clOptions.setIrange(0);
            else
                clOptions.setIrange(Double.parseDouble(clo.getOptionValue("ir").trim()));
        }
        catch (Exception e) {
            System.out.println("Invalid command line options found. Please retry.\n");
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void parseTab() {
        String out = new String();
        Hit curHit = null;
        Query curQuery = null;
        BufferedReader br = null;
    
        try {
            String line;
            br = new BufferedReader(new FileReader(clOptions.getFilename()));
            while ((line = br.readLine()) != null) {
                String[] row = line.split("/t");
                if(curQuery == null)
                    curQuery = new Query();
                else if(curQuery.getID() != row[0]) {
                    if(curQuery.getHits().size() != 0) {
                        curQuery.setHits(orderHits(curQuery));
                        for(Hit hit : curQuery.getHits()) {
                            out = curQuery.getID() + "\t" + hit.getId() + "\t" + hit.getIdentity() + "\t" + hit.getAlign_len() +
                            "\t" + hit.getMismatch() + "\t" + hit.getGapopen() + "\t" + hit.getQuery_start() + "\t" + hit.getQuery_end() + 
                            "\t" + hit.getHit_start() + "\t" + hit.getHit_end() + "\t" + hit.getEvalue() + "\t" + hit.getBitscore() + "\t" + hit.getDef() + "\n";
                            hits.write(out);
                            out = curQuery.getID() + "\t" + curHit.getId() + "\n";
                            header.write(out);
                        }
                    }
                    else {
                        out = curQuery.getID() + "\tNo hits found.\n";
                        nohits.write(out);
                    }
                    curQuery = new Query();
                }

                curHit = new Hit();
                curQuery.setID(row[0]);
                curHit.setId(row[1]);
                curHit.setIdentity(row[2]);
                curHit.setAlign_len(row[3]);
                curHit.setMismatch(row[4]);
                curHit.setGapopen(row[5]);
                curHit.setQuery_start(row[6]);
                curHit.setQuery_end(row[7]);
                curHit.setHit_start(row[8]);
                curHit.setHit_end(row[9]);
                curHit.setEvalue(row[10]);
                curHit.setBitscore(row[11]); 

                if(row.length >= 13) {
                    curHit.setDef(row[12]);
                    if(clOptions.getType() == 'n')
                        curHit.setDeflevel(Integer.toString(curHit.getDef().length() - curHit.getDef().replace(";", "").length()));
                    else if(clOptions.getType() == 'p')
                        curHit.setDeflevel(Integer.toString(curHit.getDef().length() - curHit.getDef().replace(">", "").length()));
                }

                curHit.setP_identity(String.format("%.1f", 100*(Double.parseDouble(curHit.getIdentity()) / Double.parseDouble(curHit.getAlign_len()))));

                if(Double.parseDouble(curHit.getEvalue()) <= clOptions.getEvalue()
                && Double.parseDouble(curHit.getBitscore()) >= clOptions.getBitscore()
                && Integer.parseInt(curHit.getDeflevel()) >= clOptions.getDefinition()
                && Double.parseDouble(curHit.getP_identity()) >= clOptions.getIdentity())
                    curQuery.appendHit(curHit);
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Input file could not be parsed. Check the BLAST results file.\n");
            System.exit(1);
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private void parseXML() {
        try {
            File infile = new File(clOptions.getFilename());
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(infile);
            doc.getDocumentElement().normalize();
           
            NodeList n = doc.getElementsByTagName("BlastOutput");
            NodeList List = ((Element) n.item(0)).getElementsByTagName("BlastOutput_iterations");
            NodeList nList = ((Element) List.item(0)).getElementsByTagName("Iteration");
            
            for(int idx=0; idx < nList.getLength(); idx++) {
                Node nNode = nList.item(idx);
                if(nNode.getNodeType() == Node.ELEMENT_NODE) {
                    Element eElement = (Element) nNode;
                    Query curQuery = new Query();
                    curQuery.setNum(eElement.getElementsByTagName("Iteration_iter-num").item(0).getTextContent());
                    curQuery.setDef(eElement.getElementsByTagName("Iteration_query-def").item(0).getTextContent());
                    curQuery.setLength(eElement.getElementsByTagName("Iteration_query-len").item(0).getTextContent());
                    
                    NodeList l = eElement.getElementsByTagName("Iteration_hits");
                    NodeList hitList = ((Element) l.item(0)).getElementsByTagName("Hit");
                    
                    for(int idx1=0; idx1 < hitList.getLength(); idx1++) {
                        Node hitNode = hitList.item(idx1);
                        if(hitNode.getNodeType() == Node.ELEMENT_NODE) {
                            Element hitElement = (Element) hitNode;
                            Hit curHit = new Hit();
                            curHit.setId(hitElement.getElementsByTagName("Hit_id").item(0).getTextContent());
                            curHit.setDef(hitElement.getElementsByTagName("Hit_def").item(0).getTextContent());
                            if(clOptions.getType() == 'n')
                                curHit.setDeflevel(Integer.toString(curHit.getDef().length() - curHit.getDef().replace(";", "").length()));
                            else if(clOptions.getType() == 'p')
                                curHit.setDeflevel(Integer.toString(curHit.getDef().length() - curHit.getDef().replace(">", "").length()));
                            curHit.setAccession(hitElement.getElementsByTagName("Hit_accession").item(0).getTextContent());
                            curHit.setLength(hitElement.getElementsByTagName("Hit_len").item(0).getTextContent());
                            
                            int count = 1;
                            NodeList h = hitElement.getElementsByTagName("Hit_hsps");
                            NodeList hspList = ((Element) h.item(0)).getElementsByTagName("Hsp");
                            
                            for(int idx2=0; idx2 < hspList.getLength(); idx2++) {
                                Node hspNode = hspList.item(idx2);
                                if(hspNode.getNodeType() == Node.ELEMENT_NODE) {
                                    Element hspElement = (Element) hspNode;
                                    
                                    if(count > 1){
                                        Hit newHit = new Hit();
                                        newHit.setId(curHit.getId());
                                        newHit.setDef(curHit.getDef());
                                        newHit.setDeflevel(curHit.getDeflevel());
                                        newHit.setAccession(curHit.getAccession());
                                        newHit.setLength(curHit.getLength());
                                        curHit = newHit;
                                    }
                                    
                                    curHit.setBitscore(hspElement.getElementsByTagName("Hsp_bit-score").item(0).getTextContent());
                                    curHit.setScore(hspElement.getElementsByTagName("Hsp_score").item(0).getTextContent());
                                    curHit.setEvalue(hspElement.getElementsByTagName("Hsp_evalue").item(0).getTextContent());
                                    curHit.setQuery_start(hspElement.getElementsByTagName("Hsp_query-from").item(0).getTextContent());
                                    curHit.setQuery_end(hspElement.getElementsByTagName("Hsp_query-to").item(0).getTextContent());
                                    curHit.setHit_start(hspElement.getElementsByTagName("Hsp_hit-from").item(0).getTextContent());
                                    curHit.setHit_end(hspElement.getElementsByTagName("Hsp_hit-to").item(0).getTextContent());
                                    curHit.setQuery_frame(hspElement.getElementsByTagName("Hsp_query-frame").item(0).getTextContent());
                                    curHit.setIdentity(hspElement.getElementsByTagName("Hsp_identity").item(0).getTextContent());
                                    curHit.setAlign_len(hspElement.getElementsByTagName("Hsp_align-len").item(0).getTextContent());
                                    curHit.setPositive(hspElement.getElementsByTagName("Hsp_positive").item(0).getTextContent());
                                    count++;

                                    curHit.setP_identity(String.format("%.1f", 100*(Double.parseDouble(curHit.getIdentity()) / Double.parseDouble(curHit.getAlign_len()))));
                                    curHit.setP_conserved(String.format("%.1f", 100*(Double.parseDouble(curHit.getPositive()) / Double.parseDouble(curHit.getAlign_len()))));
                                    
                                    if(Double.parseDouble(curHit.getEvalue()) <= clOptions.getEvalue()
                                    && Double.parseDouble(curHit.getBitscore()) >= clOptions.getBitscore()
                                    && Integer.parseInt(curHit.getDeflevel()) >= clOptions.getDefinition()
                                    && Double.parseDouble(curHit.getP_identity()) >= clOptions.getIdentity())
                                        curQuery.appendHit(curHit);
                                }
                            }
                        }
                    }

                    String out = new String();
                    if(curQuery.getHits().size() != 0) {
                        curQuery.setHits(orderHits(curQuery));
                        for(Hit hit : curQuery.getHits()) {
                            out = curQuery.getDef() + "\t" + curQuery.getLength() + "\t" + hit.getAccession() + "\t" + hit.getLength() +
                            "\t" + hit.getDef() + "\t" + hit.getEvalue() + "\t" + hit.getBitscore() + "\t" + hit.getQuery_frame() + 
                            "\t" + hit.getQuery_start() + "\t" + hit.getQuery_end() + "\t" + hit.getHit_start() + "\t" + hit.getHit_end() + 
                            "\t" + hit.getP_conserved() + '\t' + hit.getP_identity() + "\n";
                            hits.write(out);
                            out = curQuery.getDef() + "\t" + hit.getDef() + "\n";
                            nohits.write(out);
                        }
                    }
                    else {
                        out = curQuery.getDef() + "\tNo hits found.\n";
                        header.write(out);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("XML file could not be parsed. Check the BLAST results file.\n");
            System.exit(1);
        }
    }

    private ArrayList<Hit> orderHits(Query curQuery) {
        ArrayList<Hit> 
            hitList = curQuery.getHits(), 
            top_hits = new ArrayList<Hit>()
        ;
    
        if(clOptions.getOrder() == 'e')
            hitList.sort((h1, h2) -> { return (((Double) Double.parseDouble(h1.getEvalue())).compareTo((Double) Double.parseDouble(h2.getEvalue()))); });
        else if(clOptions.getOrder() == 'b')
            hitList.sort((h1, h2) -> { return -(((Double) Double.parseDouble(h1.getBitscore())).compareTo((Double) Double.parseDouble(h2.getBitscore()))); });
        else if(clOptions.getOrder() == 'i')
            hitList.sort((h1, h2) -> { return -(((Double) Double.parseDouble(h1.getIdentity())).compareTo((Double) Double.parseDouble(h2.getIdentity()))); });
        else if(clOptions.getOrder() == 'd')
            hitList.sort((h1, h2) -> { return -(((Integer) Integer.parseInt(h1.getDeflevel())).compareTo((Integer) Integer.parseInt(h2.getDeflevel()))); });

        if(clOptions.getErange() > 0) {
            double acceptval = Double.parseDouble(hitList.get(0).getEvalue()) + clOptions.getErange();
            for(Hit hit : hitList) {
                if(Double.parseDouble(hit.getEvalue()) <= acceptval)
                    top_hits.add(hit);
                else 
                    break;
            }
            top_hits.sort((h1, h2) -> { return -(((Integer) Integer.parseInt(h1.getDeflevel())).compareTo((Integer) Integer.parseInt(h2.getDeflevel()))); });
        }

        else if(clOptions.getBrange() > 0) {
            double acceptval = Double.parseDouble(hitList.get(0).getBitscore()) - clOptions.getBrange();
            for(Hit hit : hitList) {
                if(Double.parseDouble(hit.getBitscore()) >= acceptval)
                    top_hits.add(hit);
                else 
                    break;
            }
            top_hits.sort((h1, h2) -> { return -(((Integer) Integer.parseInt(h1.getDeflevel())).compareTo((Integer) Integer.parseInt(h2.getDeflevel()))); });
        }

        else if(clOptions.getIrange() > 0) {
            double acceptval = Double.parseDouble(hitList.get(0).getIdentity()) - clOptions.getIrange();
            for(Hit hit : hitList) {
                if(Double.parseDouble(hit.getEvalue()) >= acceptval)
                    top_hits.add(hit);
                else 
                    break;
            }
            top_hits.sort((h1, h2) -> { return -(((Integer) Integer.parseInt(h1.getDeflevel())).compareTo((Integer) Integer.parseInt(h2.getDeflevel()))); });
        }

        else 
            top_hits = hitList;

        if(clOptions.getNumber() > 0) {
            hitList.subList(clOptions.getNumber()-1, hitList.size()-1).clear();
        }

        return hitList;
    }
}