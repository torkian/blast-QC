package BLAST_QC_JAVA;

import java.util.ArrayList;

public class Query {

    private String ID = null;
    private String num = null;
    private String def = null;
    private String length = null;
    private ArrayList<Hit> hits = new ArrayList<Hit>();


    public String getID() {
        return ID;
    }

    public void setID(String iD) {
        ID = iD;
    }

    public String getNum() {
        return this.num;
    }

    public void setNum(String num) {
        this.num = num;
    }

    public String getDef() {
        return this.def;
    }

    public void setDef(String def) {
        this.def = def;
    }

    public String getLength() {
        return this.length;
    }

    public void setLength(String length) {
        this.length = length;
    }

    public ArrayList<Hit> getHits() {
        return this.hits;
    }

    public void setHits(ArrayList<Hit> hits) {
        this.hits = hits;
    }

    public void appendHit(Hit h) {
        this.hits.add(h);
    }

}