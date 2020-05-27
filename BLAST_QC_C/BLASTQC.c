#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#include "libxml/parser.h"
#include "BLASTQC.h"
#include "vector.h"

char CMP_BY = 'e'; 

typedef struct Hit {
    char* num;
    char* id;
    char* def;
    int deflevel;
    char* accession;
    char* length;
    char* mismatch;
    char* gapopen;
    double bitscore;
    char* score;
    double evalue;
    char* query_start;
    char* query_end;
    char* hit_start;
    char* hit_end;
    char* query_frame;
    int positive;
    int identity;
    int align_len;
    double p_identity;
    double p_conserved;
} Hit;

typedef struct Query {
    char* id;
    char* num;
    char* def;
    char* length;
    vector hitlist;
} Query;

bool validateArgs();
int cmpfunc (const Hit* a, const Hit* b);
int partition_r(vector* hitlist, int low, int high);
int partition (vector* hitlist, int low, int high);
void quickSort(vector* hitlist, int low, int high);
void orderHits(vector* hitlist);
const char* getFeild(char* line, int idx);
void parseTab(FILE* hits, FILE* nohits, FILE* header);
void parseHsp(xmlDocPtr doc, xmlNodePtr cur, Query* curQuery, Hit* curHit);
void parseHit(xmlDocPtr doc, xmlNodePtr cur, Query* curQuery);
void parseQuery(xmlDocPtr doc, xmlNodePtr cur, FILE* hits, FILE* nohits, FILE* header);
void parseXML(FILE* hits, FILE* nohits, FILE* header);

int main(int argc, char **argv) {
    if(!validateArgs())
        return 1;

    srand(time(NULL));
    FILE* hits = fopen(strcat(strdup(OUTFILE), ".hits.txt"), "w");
    FILE* nohits = fopen(strcat(strdup(OUTFILE), ".nohits.txt"), "w");
    FILE* header = fopen(strcat(strdup(OUTFILE), ".hits.header"), "w");

    if(strcmp(FILEFORMAT, "XML") == 0) {
        fprintf(
            hits, "query_name\tquery_length\taccession_number\tsubject_length\t"
            "subject_description\tE_value\tbit_score\tframe\tquery_start\tquery_end\t"
            "hit_start\thit_end\t%%conserved\t%%identity\n"
        );
        fprintf(nohits, "query_name\n");
        fprintf(header, "query_name\tsubject_description\n");
        parseXML(hits, nohits, header);
    }
           
    else if(strcmp(FILEFORMAT, "tab") == 0) {
        fprintf(
            hits, "qseqid\tsseqid\tpident\tlength\tmismatch\t"
            "gapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
        );
        fprintf(nohits, "qseqid\n");
        fprintf(header, "qseqid\tsseqid\n");
        parseTab(hits, nohits, header);
    }

    fclose(hits);
    fclose(nohits);
    fclose(header);

    return 0;
}


bool validateArgs() {
    bool hasReqArgs = true;

    if(strlen(FILENAME) < 1) {
        printf("Must specify the BLAST results infile path.");
        hasReqArgs = false;
    }

    if(strlen(FILEFORMAT) < 1) {
        printf("Must specify a result file format (XML or tab).\n");
        hasReqArgs = false;
    }

    if(TYPE != 'n' && TYPE != 'p') {
        printf("Must specify the type of BLAST you are using (n or p).\n");
        hasReqArgs = false;
    }

    if(ERANGE > 0 && ORDER != 'e') {
        printf("erange cannot be used. Must order by evalue if this functionality is desired.");
        hasReqArgs = false;
    }

    if(BRANGE > 0 && ORDER != 'b') {
        printf("brange cannot be used. Must order by bitscore if this functionality is desired.");
        hasReqArgs = false;
    }

    if(IRANGE > 0 && ORDER != 'i') {
        printf("irange cannot be used. Must order by identity if this functionality is desired.");
        hasReqArgs = false;
    }

    if(!hasReqArgs) {
        printf("Invalid Args. Please check the BLASTQC.h header file and retry.\n");
    }

    if(ORDER != '\0')
        CMP_BY = ORDER;

    if(strlen(OUTFILE) < 1)
        strcpy(OUTFILE, FILENAME);
    
    return hasReqArgs;
}

//Compare 2 Hits by the input comparison value
int cmpfunc (const Hit* a, const Hit* b) {
    if(CMP_BY == 'e')
        return (a->evalue < b->evalue);
    else if(CMP_BY == 'b')
        return (a->bitscore > b->bitscore);
    else if(CMP_BY == 'i')
        return (a->p_identity > b->p_identity);
    return -1;
} 
  
/* This function takes a random element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition (vector* hitlist, int low, int high) {
    Hit* pivot = ((Hit*) vector_get(hitlist, high));    // pivot 

    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (cmpfunc(((Hit*) vector_get(hitlist, j)), pivot))
        { 
            i++;    // increment index of smaller element 
            Hit* x = ((Hit*) vector_get(hitlist, i));
            Hit* y = ((Hit*) vector_get(hitlist, j));
            vector_set(hitlist, j, x);
            vector_set(hitlist, i, y); 
        } 
    } 
    Hit* a = ((Hit*) vector_get(hitlist, i + 1));
    Hit* b = ((Hit*) vector_get(hitlist, high));
    vector_set(hitlist, high, a);
    vector_set(hitlist, i + 1, b);
 
    return (i + 1); 
} 

// A utility function to randomize partition pivot element
int partition_r(vector* hitlist, int low, int high) {  
    int pidx = low + rand()%(high - low + 1);
    Hit* a = ((Hit*) vector_get(hitlist, pidx));
    Hit* b = ((Hit*) vector_get(hitlist, high));
    vector_set(hitlist, high, a);
    vector_set(hitlist, pidx, b);

    return partition(hitlist, low, high);
}
  
// The main function that implements QuickSort 
void quickSort(vector* hitlist, int low, int high) { 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition_r(hitlist, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(hitlist, low, pi - 1); 
        quickSort(hitlist, pi + 1, high); 
    } 
} 

void orderHits(vector* hitlist) {
    quickSort(hitlist, 0, vector_total(hitlist)-1);

    vector* tophits = malloc(sizeof(vector));
    vector_init(tophits);

    if(ERANGE > 0) {
        double acceptVal = ((Hit*) vector_get(hitlist, 0))->evalue + ERANGE; 
        for(int i=0; i<vector_total(hitlist); i++) {
            if(((Hit*) vector_get(hitlist, i))->evalue <= acceptVal)
                vector_add(tophits, ((Hit*) vector_get(hitlist, i)));
            else
                break;    
        }
        CMP_BY = 'd';
        quickSort(tophits, 0, vector_total(tophits)-1);
    }

    else if(BRANGE > 0) {
        double acceptVal = ((Hit*) vector_get(hitlist, 0))->bitscore - BRANGE; 
        for(int i=0; i<vector_total(hitlist); i++) {
            if(((Hit*) vector_get(hitlist, i))->bitscore >= acceptVal)
                vector_add(tophits, ((Hit*) vector_get(hitlist, i)));
            else
                break;    
        }
        CMP_BY = 'd';
        quickSort(tophits, 0, vector_total(tophits)-1);
    }

    else if(IRANGE > 0) {
        double acceptVal = ((Hit*) vector_get(hitlist, 0))->p_identity - IRANGE; 
        for(int i=0; i<vector_total(hitlist); i++) {
            if(((Hit*) vector_get(hitlist, i))->p_identity >= acceptVal)
                vector_add(tophits, ((Hit*) vector_get(hitlist, i)));
            else
                break;    
        }
        CMP_BY = 'd';
        quickSort(tophits, 0, vector_total(tophits)-1);
    }

    else 
        tophits = hitlist;


    vector_free(hitlist);
    hitlist = tophits;
}

const char* getFeild(char* line, int idx) {
    const char* tok;
    for (tok = strtok(line, "\t");
            tok && *tok;
            tok = strtok(NULL, "\t"))
    {
        if (!idx)
            return tok;
    }
    return NULL;
}

void parseTab(FILE* hits, FILE* nohits, FILE* header) {
    FILE* infile = fopen(FILENAME, "r");
    char* row;

    Query* curQuery = NULL;
    Hit* curHit;

    while(fgets(row, BUFFER_SIZE, infile)) {
        if(curQuery == NULL) {
            curQuery = malloc(sizeof(Query));
            vector_init(&curQuery->hitlist);
        }

        else if(strcmp(curQuery->id, getFeild(strdup(row), 0)) != 0) {
            if(vector_total(&curQuery->hitlist)) {
                orderHits(&curQuery->hitlist);
                for(int i=0; i<vector_total(&curQuery->hitlist); i++) {
                    fprintf(
                        hits, "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%s\n", 
                        curQuery->id, ((Hit*) vector_get(&curQuery->hitlist, i))->id, ((Hit*) vector_get(&curQuery->hitlist, i))->identity,
                        ((Hit*) vector_get(&curQuery->hitlist, i))->align_len, ((Hit*) vector_get(&curQuery->hitlist, i))->mismatch, ((Hit*) vector_get(&curQuery->hitlist, i))->gapopen,
                        ((Hit*) vector_get(&curQuery->hitlist, i))->query_start, ((Hit*) vector_get(&curQuery->hitlist, i))->query_end, ((Hit*) vector_get(&curQuery->hitlist, i))->hit_start,
                        ((Hit*) vector_get(&curQuery->hitlist, i))->hit_end, ((Hit*) vector_get(&curQuery->hitlist, i))->evalue, ((Hit*) vector_get(&curQuery->hitlist, i))->bitscore,
                        ((Hit*) vector_get(&curQuery->hitlist, i))->def
                    );
                    fprintf(header, "%s\t%s\n", curQuery->id, ((Hit*) vector_get(&curQuery->hitlist, i))->id);
                }
            }
            else 
                fprintf(nohits, "%s\tNo hits found.\n", curQuery->id);
            
            vector_free(&curQuery->hitlist);
            free(curQuery);
            curQuery = malloc(sizeof(Query));
            vector_init(&curQuery->hitlist);
        }

        curHit = malloc(sizeof(Hit));
        strcpy(curQuery->id, getFeild(strdup(row), 0));
        strcpy(curHit->id, getFeild(strdup(row), 1));
        sscanf(getFeild(strdup(row), 2), "%d", &curHit->identity);
        sscanf(getFeild(strdup(row), 3), "%d", &curHit->align_len);
        strcpy(curHit->mismatch, getFeild(strdup(row), 4));
        strcpy(curHit->gapopen, getFeild(strdup(row), 5));
        strcpy(curHit->query_start, getFeild(strdup(row), 6));
        strcpy(curHit->query_end, getFeild(strdup(row), 7));
        strcpy(curHit->hit_start, getFeild(strdup(row), 8));
        strcpy(curHit->hit_end, getFeild(strdup(row), 9));
        sscanf(getFeild(strdup(row), 10), "%lf", &curHit->evalue);
        sscanf(getFeild(strdup(row), 11), "%lf", &curHit->bitscore);
        char str[32];
        sprintf(str, "%.1lf", 100*(curHit->identity/(double)curHit->align_len));
        sscanf(str, "%lf", &curHit->p_identity);

        vector_add(&curQuery->hitlist, curHit);
    }

    free(curQuery);
}

void parseHsp(xmlDocPtr doc, xmlNodePtr cur, Query* curQuery, Hit* curHit) {
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
        if(cur->type == XML_ELEMENT_NODE) {
            if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_bit-score"))) {
                sscanf((char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1), "%lf", &curHit->bitscore);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_evalue"))) {
                sscanf((char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1), "%lf", &curHit->evalue);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_query-from"))) {
                curHit->query_start = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_query-to"))) {
                curHit->query_end = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_hit-from"))) {
                curHit->hit_start = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_hit-to"))) {
                curHit->hit_end = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_query-frame"))) {
                curHit->query_frame = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_identity"))) {
                sscanf((char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1), "%d", &curHit->identity);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_positive"))) {
                sscanf((char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1), "%d", &curHit->positive);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp_align-len"))) {
                sscanf((char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1), "%d", &curHit->align_len);
            }
        }
        cur = cur->next;
    }

    char str[32];
    sprintf(str, "%.1lf", 100*(curHit->identity/(double)curHit->align_len));
    sscanf(str, "%lf", &curHit->p_identity);
    sprintf(str, "%.1lf", 100*(curHit->positive/(double)curHit->align_len));
    sscanf(str, "%lf", &curHit->p_conserved);
}

void parseHit(xmlDocPtr doc, xmlNodePtr cur, Query* curQuery) {
    int count = 0;
    Hit* curHit = malloc(sizeof(Hit));
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
        if(cur->type == XML_ELEMENT_NODE) {
            if((!xmlStrcmp(cur->name, (const xmlChar *) "Hit_id"))) {
                curHit->id = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hit_def"))) {
                curHit->def = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hit_accession"))) {
                curHit->accession = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hit_len"))) {
                curHit->length = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hit_hsps"))) {
                cur = cur->xmlChildrenNode;
                continue;
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hsp"))) {
                count++;
                if(count > 1) {
                    Hit* newHit = malloc(sizeof(Hit));
                    newHit->deflevel = curHit->deflevel;
                    strcpy(newHit->id, curHit->id);
                    strcpy(newHit->def, curHit->def);
                    strcpy(newHit->accession, curHit->accession);
                    strcpy(newHit->length, curHit->length);
                    curHit = newHit;
                }
                parseHsp(doc, cur, curQuery, curHit);
                vector_add(&curQuery->hitlist, curHit);
            }
        }
        cur = cur->next;
    }
    
}

void parseQuery(xmlDocPtr doc, xmlNodePtr cur, FILE* hits, FILE* nohits, FILE* header) {
    Query* curQuery = malloc(sizeof(Query));
    vector_init(&curQuery->hitlist);
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
        if(cur->type == XML_ELEMENT_NODE) {
            if((!xmlStrcmp(cur->name, (const xmlChar *) "Iteration_iter-num"))) {
                curQuery->num = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Iteration_query-def"))) {
                curQuery->def = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Iteration_query-len"))) {
                curQuery->length = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Iteration_hits"))) {
                cur = cur->xmlChildrenNode;
                continue;
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Hit"))) {
                parseHit(doc, cur, curQuery);
            }
        }
        cur = cur->next;
    }

    if(vector_total(&curQuery->hitlist)) {
        orderHits(&curQuery->hitlist);
        for(int i=0; i<vector_total(&curQuery->hitlist); i++) {
            fprintf(
                hits, "%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\n", 
                curQuery->def, curQuery->length, ((Hit*) vector_get(&curQuery->hitlist, i))->accession, ((Hit*) vector_get(&curQuery->hitlist, i))->length,
                ((Hit*) vector_get(&curQuery->hitlist, i))->def, ((Hit*) vector_get(&curQuery->hitlist, i))->evalue, ((Hit*) vector_get(&curQuery->hitlist, i))->bitscore,
                ((Hit*) vector_get(&curQuery->hitlist, i))->query_frame, ((Hit*) vector_get(&curQuery->hitlist, i))->query_start, ((Hit*) vector_get(&curQuery->hitlist, i))->query_end,
                ((Hit*) vector_get(&curQuery->hitlist, i))->hit_start, ((Hit*) vector_get(&curQuery->hitlist, i))->hit_end, ((Hit*) vector_get(&curQuery->hitlist, i))->p_conserved, 
                ((Hit*) vector_get(&curQuery->hitlist, i))->p_identity
            );
            fprintf(header, "%s\t%s\n", curQuery->def, ((Hit*) vector_get(&curQuery->hitlist, i))->def);
        }
    }
    else 
        fprintf(nohits, "%s\tNo hits found.\n", curQuery->def);

    free(curQuery);
}

void parseXML(FILE* hits, FILE* nohits, FILE* header) {
    xmlDocPtr doc = xmlParseFile(FILENAME);

    if(doc == NULL) {
        printf("Failed to parse the BLAST XML results file %s. Please try again.\n", FILENAME);
        exit(-1);
    }

    xmlNodePtr cur = xmlDocGetRootElement(doc);

    if(cur == NULL) {
        printf("BLAST Results XML file: %s is empty. Please try again.\n", FILENAME);
        exit(1);
    }

    Query* curQuery = NULL;
    Hit* curHit = NULL;

    if(xmlStrcmp (cur->name, (const xmlChar *) "BlastOutput")){
        printf ("XML document: %s is not a BLAST Output file. Root node != BlastOutput", FILENAME);
        exit(2);
    }

    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
        if(cur->type == XML_ELEMENT_NODE) {
            if((!xmlStrcmp(cur->name, (const xmlChar *) "BlastOutput_iterations"))) {
                cur = cur->xmlChildrenNode;
                continue;
            }
            else if((!xmlStrcmp(cur->name, (const xmlChar *) "Iteration"))) {
                parseQuery(doc, cur, hits, nohits, header);
            }
        }
        cur = cur->next;
    }

    xmlFree(doc);
}
