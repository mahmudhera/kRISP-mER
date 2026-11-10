//
//  main.cpp
//  kmer
//

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>
using namespace std;

#include <pthread.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "kmer.h"



#define MAX_REC_LEN 10240

int noCases, noControls;

int kmerLength=23;

int NUM_THREADS=2;

pthread_mutex_t readFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t caseOutFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t controlOutFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t refFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t refMultFile_mutex = PTHREAD_MUTEX_INITIALIZER;


pthread_mutex_t hashTable_mutex = PTHREAD_MUTEX_INITIALIZER;

#pragma pack(push, 1)
class Kmer
{
public:
	long long int kmer;
	unsigned long caseCounts;
    	unsigned long controlCounts;
	
    	Kmer(int noCases, int noControls);
    	~Kmer();
    	void show();
	void freeMemory();
};
#pragma pack(pop)

class KeyVal
{
public:
	long long int val;
	long count;
};

class HashTable
{
public:
	long long int totalKmers;
    	long long int totalTests;
    	vector <Kmer *> kmers[HASH_TABLE_LENGTH];
	pthread_mutex_t hashTableBuckets_mutex[HASH_TABLE_LENGTH];

	long long int totalKmers_Bucket[HASH_TABLE_LENGTH];	

	
	long long int totalKmerCountsCase;	
	long long int totalKmerCountsControl;
    	HashTable();
    	void insertKmer(long long int val, int count, int isCase);
    	void show();
    	void dumpKmers();
    
};

HashTable *ht;


long long int getInt(char *s)
{
	long long int val=0;
	int i=0;
	char ch;
	while(s[i])
	{
		val=val<<2;
		
		ch=s[i];
		if(ch=='A')
		{
			val=val|0;
		}
		else if(ch=='C')
		{
			val=val|1;
		}
		else if(ch=='G')
		{
			val=val|2;
		}
		else
		{
			val=val|3;
		}
		i++;
        
	}
	return val;
}

char * getKmer(long long int val, char *kmer, int kmerLength)
{
    
	int temp=0;
	for(int i=kmerLength-1;i>=0;i--)
	{
		temp=val&3;
		val=val>>2;
		if(temp==0)
			kmer[i]='A';
		else if(temp==1)
			kmer[i]='C';
		else if(temp==2)
			kmer[i]='G';
		else if(temp==3)
			kmer[i]='T';
	}
	kmer[kmerLength]='\0';
	return kmer;
}


unsigned long int getHash(unsigned long long int key)
{
    /*
     key = (~key) + (key << 18); // key = (key << 18) - key - 1;
     key = key ^ (key >> 31);
     key = key * 21; // key = (key + (key << 2)) + (key << 4);
     key = key ^ (key >> 11);
     key = key + (key << 6);
     key = key ^ (key >> 22);
     */
	return (unsigned long int) key;
}




Kmer::Kmer (int noCases, int noControls)
{
	caseCounts=0;
    
	
   	 controlCounts=0;
}

Kmer::~Kmer()
{
//    delete [] caseCounts;
//    delete [] controlCounts;
}

void Kmer::freeMemory()
{
//	delete [] caseCounts;
//	delete [] controlCounts;
//	free(caseCounts);
//	free(controlCounts);
}


void Kmer::show()
{
    char kmerString[100];
    
    cout<<getKmer(kmer,kmerString,23)<<" ";
    
    cout<<caseCounts<<" ";
        
    cout<<controlCounts<<endl;
        
    
}


HashTable::HashTable()
{
    totalKmerCountsCase=0;

    totalKmerCountsControl=0;

	for(int i=0;i<HASH_TABLE_LENGTH;i++)
	{
		pthread_mutex_init(&hashTableBuckets_mutex[i],NULL);
		totalKmers_Bucket[i]=0;
	}
    
    totalKmers=0;
    totalTests=0;
}


void HashTable::insertKmer(long long int val, int count, int isCase)
{
	unsigned long int index=getHash(val) % HASH_TABLE_LENGTH;
	int found=0;

	pthread_mutex_lock(&hashTableBuckets_mutex[index]);
    
	for(int i=0;i<kmers[index].size();i++)
	{
		if(val==kmers[index][i]->kmer)
        	{

            		if(isCase==1)
            		{
                		kmers[index][i]->caseCounts=count;
            		}
            		else
            		{
                		kmers[index][i]->controlCounts=count;
            		}
			found=1;
			break;
		}
	}
	if(found==0)
	{
		Kmer* km=new Kmer(noCases,noControls);
		km->kmer=val;
		
        	if(isCase==1)
        	{
            		km->caseCounts=count;
        	}
        	else
        	{
            		km->controlCounts=count;
        	}
        
		kmers[index].push_back(km);
       
		totalKmers_Bucket[index]++;
        	

	}
	pthread_mutex_unlock(&hashTableBuckets_mutex[index]);
 
}




FILE *caseFile;
FILE *controlFile;
FILE *refFile;
FILE *refMultFile;


void * dump_thread(void *threadid)
{
	long tid;
   	tid = (long)threadid;

	char kmerString[100];
    
    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {
            if(ht->kmers[i][j]->caseCounts>0 && ht->kmers[i][j]->controlCounts==0)
            {
          

			pthread_mutex_lock(&caseOutFile_mutex);

                        fprintf(caseFile,"%s\t%d\t%d\n",getKmer(ht->kmers[i][j]->kmer, kmerString, 23),(int)ht->kmers[i][j]->caseCounts,(int)ht->kmers[i][j]->controlCounts);

			pthread_mutex_unlock(&caseOutFile_mutex);


             
                             
            }
	    if(ht->kmers[i][j]->caseCounts>1 && ht->kmers[i][j]->controlCounts==1)
            {
          

			pthread_mutex_lock(&controlOutFile_mutex);

                        fprintf(controlFile,"%s\t%d\t%d\n",getKmer(ht->kmers[i][j]->kmer, kmerString, 23),(int)ht->kmers[i][j]->caseCounts,(int)ht->kmers[i][j]->controlCounts);

			pthread_mutex_unlock(&controlOutFile_mutex);


             
                             
            }

	    if(ht->kmers[i][j]->caseCounts==0 && ht->kmers[i][j]->controlCounts>0)
            {
          

			pthread_mutex_lock(&refFile_mutex);

                        fprintf(refFile,"%s\t%d\t%d\n",getKmer(ht->kmers[i][j]->kmer, kmerString, 23),(int)ht->kmers[i][j]->caseCounts,(int)ht->kmers[i][j]->controlCounts);

			pthread_mutex_unlock(&refFile_mutex);


             
                             
            }
	    if(ht->kmers[i][j]->caseCounts==1 && ht->kmers[i][j]->controlCounts>1)
            {
          

			pthread_mutex_lock(&refMultFile_mutex);

                        fprintf(refMultFile,"%s\t%d\t%d\n",getKmer(ht->kmers[i][j]->kmer, kmerString, 23),(int)ht->kmers[i][j]->caseCounts,(int)ht->kmers[i][j]->controlCounts);

			pthread_mutex_unlock(&refMultFile_mutex);


             
                             
            }
		
	
        }
        
    }
    
    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {
		ht->kmers[i][j]->freeMemory();
		delete ht->kmers[i][j];
        }
        ht->kmers[i].clear();
    }


	pthread_exit(NULL);

}


void HashTable::dumpKmers()
{
    caseFile=fopen("present_v_absent.txt","a");
    controlFile=fopen("multiple_v_single.txt","a");
	
	refFile=fopen("ref_present_v_absent.txt","a");
	refMultFile=fopen("ref_multiple_v_single.txt","a");

	
	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for(t=0; t<NUM_THREADS; t++)
	{
		  rc = pthread_create(&threads[t], NULL, dump_thread, (void *)t);
		  if (rc){
			 exit(-1);
		  }
	}
	
	for(t=0; t<NUM_THREADS; t++) 
	{
		rc = pthread_join(threads[t], &status);
		if (rc) 
	  	{
         		exit(-1);
      		}
      	}

	fclose(caseFile);
	fclose(controlFile);
	fclose(refFile);
	fclose(refMultFile);


    
}

void HashTable::show()
{
/*    
    for(int i=0;i<HASH_TABLE_LENGTH;i++)
    {
        for(int j=0;j<kmers[i].size();j++)
        {
            kmers[i][j]->show();
        }
        
    }
*/
    cout<<totalKmers<<endl;
    cout<<totalTests<<endl;
}

		

void getKeyVal(char *s, KeyVal* kv)
{
	long long int	val=0;
	int countVal=0;
	int i=0;
	char ch;
	while(i<KMER_LENGTH)
	{
		
		ch=s[i++];
		
		val=val<<2|bases[ch];		
        
	}
	i++;
	kv->val=val;
	while(1)
	{
		ch=s[i];
		if(ch=='\0'||ch=='\n')
		{	
			break;
		}
		countVal=countVal*10+ch-'0';
		i++;
	        
	}
	kv->count=countVal;

}

FILE * kmerFilesCases;
FILE * kmerFilesControls;	
long long int valsCases;
long long int valsControls;
unsigned long countsCases;
unsigned long countsControls;

struct ThreadArg
{
	long long int valBar;
	int threadID;
};

void * readCases(void *threadid)
{
	ThreadArg *ta=(ThreadArg *)threadid;
	long long int valBar=ta->valBar;
   	int threadNo=ta->threadID;
	long long int val;
	int count;

	char kmerString[24];

	char *line= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

	for(int i=threadNo;i<noCases;i+=NUM_THREADS/2)
    		{
			if(valsCases!=-1 && valsCases<valBar)
			{
				ht->insertKmer(valsCases, countsCases, 1);
				valsCases=-1;
				countsCases=-1;
			}
			if(valsCases==-1)
			{
       		//	while(fscanf(kmerFilesCases,"%lld %d\n",&val,&count)!=EOF)
			while(fscanf(kmerFilesCases,"%s %d\n",kmerString,&count)!=EOF)
        		{
			
				val=getInt(kmerString);
											
				if(val<valBar)
				{
 			           	ht->insertKmer(val, count, 1);
            			}
				else
				{
					valsCases=val;
					countsCases=count;
					break;
				}
        		}
 			}
    		}


	delete []line;
	pthread_exit(NULL);

}

void * readControls(void *threadid)
{
	ThreadArg *ta=(ThreadArg *)threadid;
	long long int valBar=ta->valBar;
   	int threadNo=ta->threadID;
	long long int val;
	int count;
	
	char kmerString[24];


	char *line= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);


	for(int i=threadNo;i<noControls;i+=NUM_THREADS/2)
    		{
			if(valsControls!=-1 && valsControls<valBar)
			{
				ht->insertKmer(valsControls, countsControls, 0);
				valsControls=-1;
				countsControls=-1;
			}
			if(valsControls==-1)
			{
       	//		while(fscanf(kmerFilesControls,"%lld %d\n",&val,&count)!=EOF)
        		while(fscanf(kmerFilesControls,"%s %d\n",kmerString,&count)!=EOF)
        		{
			
				val=getInt(kmerString);
			

            			if(val<valBar)
				{
 			           	ht->insertKmer(val, count, 0);
            			}
				else
				{
					valsControls=val;
					countsControls=count;
					break;
				}
			
        		}
 			}
    		}

	delete []line;
	pthread_exit(NULL);

}



int main(int argc, const char * argv[])
{

    noCases=1;
    noControls=1;
    
 
    ht=new HashTable();

	
	FILE *caseFile=fopen("present_v_absent.txt","w");
    	FILE *controlFile=fopen("multiple_v_single.txt","w");
	FILE *refFile=fopen("ref_present_v_absent.txt","w");
	FILE *refMultFile=fopen("ref_multiple_v_single.txt","w");
	

	fclose(caseFile);
	fclose(controlFile);
	fclose(refFile);
	fclose(refMultFile);



    
	char *kmerFilename;
    	kmerFilename=new char[5000];
    

	char *line= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);



    
	    char *temp;
    	char kmerString[100];
    	unsigned long count;
	
	long long int valMax=0x3FFFFFFFFFFFF;
	long long int valBar=0;
	long long int val=0;
	long long int valInc=0x0010000000000;

	FILE *sortedFile=fopen("case_sorted_files.txt","r");
	 
	for(int i=0;i<noCases;i++)
    	{
		fscanf(sortedFile,"%s\n",kmerFilename);
       		kmerFilesCases=fopen(kmerFilename,"r");

		if(kmerFilesCases==NULL)
		{
			cout<<kmerFilename<<" file doesn't exist"<<endl;
		}

		valsCases=-1;
		countsCases=-1;
		
	}
	sortedFile=fopen("control_sorted_files.txt","r");
	for(int i=0;i<noControls;i++)
    	{
		fscanf(sortedFile,"%s\n",kmerFilename);
	       	kmerFilesControls=fopen(kmerFilename,"r");

		if(kmerFilesControls==NULL)
		{
			cout<<kmerFilename<<" file doesn't exist"<<endl;
		}

		valsControls=-1;
		countsControls=-1;
		
	}

	
	ThreadArg *thArgsCase[NUM_THREADS/2];
	ThreadArg *thArgsControl[NUM_THREADS/2];

	while(valBar<valMax)
	{
		valBar+=valInc;


		pthread_t caseThreads[NUM_THREADS/2];
		pthread_t controlThreads[NUM_THREADS/2];
		int rc;
		void *status;
		for(int i=0;i<NUM_THREADS/2;i++)
		{
			thArgsCase[i]=new ThreadArg;
			thArgsCase[i]->valBar=valBar;
			thArgsCase[i]->threadID=i;
			rc = pthread_create(&caseThreads[i], NULL, readCases, (void *)thArgsCase[i]);
			if (rc)
			{
				 exit(-1);
			}
		}
		for(int i=0;i<NUM_THREADS/2;i++)
		{
			thArgsControl[i]=new ThreadArg;
			thArgsControl[i]->valBar=valBar;
			thArgsControl[i]->threadID=i;
			rc = pthread_create(&controlThreads[i], NULL, readControls, (void *)thArgsControl[i]);
			if (rc)
			{
				 exit(-1);
			}
		}
		
		for(int i=0;i<NUM_THREADS/2;i++)
		{
			rc = pthread_join(caseThreads[i], &status);
			delete thArgsCase[i];
			if (rc) 
	  		{
         			exit(-1);
      			}
		}
      		for(int i=0;i<NUM_THREADS/2;i++)
		{
			rc = pthread_join(controlThreads[i], &status);
			delete thArgsControl[i];
			if (rc) 
	  		{
         			exit(-1);
      			}
		}


		
	        	
    		ht->dumpKmers();
		
		cout<<valBar<<endl;
 
	}

	for(int i=0;i<HASH_TABLE_LENGTH;i++)
	{
		ht->totalKmers+=ht->totalKmers_Bucket[i];
		
	}

    	ht->show();







    
    return 0;
}




