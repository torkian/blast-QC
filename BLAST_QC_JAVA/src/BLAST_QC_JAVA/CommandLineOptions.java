package BLAST_QC_JAVA;

public class CommandLineOptions {
    private String filename;
    private String fileformat;
    private String outfile;
    private char type;
    private int number;
    private double evalue;
    private double bitscore;
    private double identity;
    private int definition;
    private char order;
    private double erange;
    private double brange;
    private double irange;


    public String getFilename() {
        return this.filename;
    }

    public void setFilename(String filename) {
        this.filename = filename;
    }

    public String getFileformat() {
        return this.fileformat;
    }

    public void setFileformat(String fileformat) {
        this.fileformat = fileformat;
    }

    public String getOutfile() {
        return this.outfile;
    }

    public void setOutfile(String outfile) {
        this.outfile = outfile;
    }

    public char getType() {
        return this.type;
    }

    public void setType(char type) {
        this.type = type;
    }

    public int getNumber() {
        return this.number;
    }

    public void setNumber(int number) {
        this.number = number;
    }

    public double getEvalue() {
        return this.evalue;
    }

    public void setEvalue(double evalue) {
        this.evalue = evalue;
    }

    public double getBitscore() {
        return this.bitscore;
    }

    public void setBitscore(double bitscore) {
        this.bitscore = bitscore;
    }

    public double getIdentity() {
        return this.identity;
    }

    public void setIdentity(double identity) {
        this.identity = identity;
    }

    public int getDefinition() {
        return this.definition;
    }

    public void setDefinition(int definition) {
        this.definition = definition;
    }

    public char getOrder() {
        return this.order;
    }

    public void setOrder(char order) {
        this.order = order;
    }

    public double getErange() {
        return this.erange;
    }

    public void setErange(double erange) {
        this.erange = erange;
    }

    public double getBrange() {
        return this.brange;
    }

    public void setBrange(double brange) {
        this.brange = brange;
    }

    public double getIrange() {
        return this.irange;
    }

    public void setIrange(double irange) {
        this.irange = irange;
    }

}
