// reads_fid(fname,-1,'coilScale'); read all channels and do scaling if rawdatacorrection is true 
// reads_fid(fname,chan);   read channel chan (0 based) but no correction
// reads_fid(fname); read all channels and no correction
#include <stdlib.h> 
#include <stdio.h>
#include <string>
#include <fstream>
#include "mex.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#define MAXLINES 20000
//extern "C" int __cdecl _fseeki64(FILE *, __int64, int);
//extern "C" __int64 __cdecl _ftelli64(FILE *);

using namespace std;

void dec2bin(unsigned int d, int len, short *res)
{

	for (int i = 0; i < len; i++)
	{
		res[i] = d%2;

		d = d / 2;

	}

}

bool fexists(const char *filename) {
    
    FILE *fid=fopen(filename,"r");
    if (fid!=NULL)
    {
        fclose(fid);
        return true;
    }
     else
     return false;
 
}

long get_file_size(char* filename) // path to file
{
	
	FILE *p_file = NULL;
	p_file = fopen(filename, "rb");
	if (p_file == NULL) return 0;

	fseek(p_file, 0, SEEK_END);
	long size = ftell(p_file);
	fclose(p_file);
	return size;
	
	/*
	struct stat filestatus;
	stat(filename, &filestatus);

	return filestatus.st_size;
	*/

}

void extract_protocol(char *fname)
 {

	 string line;
	 char fpro[100];
	 sprintf(fpro, "%s.pro", fname);

	 ifstream myfile(fname);
	 if (!myfile.is_open())
	 {
		 return;
	 }

	
	 ofstream ofile(fpro);
 
 int start=0;
 char buffer[1000];
 while (getline(myfile, line))
{

	std::size_t found = line.find("ASCCONV");

	if (found != std::string::npos)
	
	{
       start=start+1;
	   printf("%s\n", buffer);
        if (start==1)
        {
           ofile << "### ASCCONV BEGIN ###\n";
            continue;
        }
        else
        {
            ofile << "### ASCCONV END ###\n";
            break;
        }
   }    
    
    
    if (start==1)
    {
		ofile << line << endl;
    }   
    
}

 ofile.close();

myfile.close();

   
 }
 
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{                 

	int chan=100;
	if (nrhs >1)
	{
		chan = (int)*mxGetPr(prhs[1]);
	}
	
	bool allchan = (nrhs == 1 || chan < 0);
    float sclr[32],scli[32];
   
    float sclft[32];
    if (nrhs>2)
    {
    
       
    char *fscale=mxArrayToString(prhs[2]);
    FILE *fid_scale=fopen(fscale,"rb");
    if (fid_scale==NULL)  
    {
        printf("File %s not found\n",fscale);
        return;
    }
     char tmp[100];
    for(int i=0;i<32;i++)
    {
     fscanf(fid_scale,"{ { %s } { %f } { %f } { %f } }\n",tmp,&sclft[i],&sclr[i],&scli[i]);
     printf("{ %s{%f}{%f}{%f} }\n",tmp,sclft[i],sclr[i],scli[i]);

    }
     fclose(fid_scale);
     
    }
    
    
char *filename=mxArrayToString(prhs[0]);

long long fsize = get_file_size(filename);


FILE *fid = fopen(filename,"rb");
int nProtHeaderLen; 
if (fid == NULL) return;

fread(&nProtHeaderLen,1,4,fid);
char buffer[100];


sprintf(buffer,"%s.header",filename);
if (fexists(buffer))
   fseek(fid, nProtHeaderLen - 4, SEEK_CUR);
else
{
    char *p=(char *)malloc(nProtHeaderLen-4);

    fread(p,1,nProtHeaderLen-4,fid);


 FILE *fid2=fopen(buffer,"w");

  for (int i=4;i<nProtHeaderLen-4-12;i++)
  {
   fprintf(fid2,"%c",p[i]);
  }
  fclose(fid2);
  free(p);

  extract_protocol(buffer);
}



bool acqEnded = false;
bool rawdataCor;
long mdhCounter = 0;
float *data; 

unsigned short *phInd, *parInd, *Rep, *slcInd, *echoInd, *chInd;

unsigned short ushSamplesInScan, ushSamplesInScanTmp;

unsigned short ushUsedChannels;
unsigned short ushUsedChannelsTmp;

long totalLine;
int oldMostSig = -1;
while (!acqEnded)
{
unsigned int ulFlagsAndDMALength;
fread(&ulFlagsAndDMALength, 1, 4, fid);

//lMeasUID = fread(fid, 1, '*int32');
fseek(fid, 16, SEEK_CUR);

unsigned int aulEvalInfoMaskMostSig;
fread(&aulEvalInfoMaskMostSig, 1,4,fid); //evaluation info mask field, first part

unsigned int aulEvalInfoMaskLeastSig;
fread(&aulEvalInfoMaskLeastSig, 1, 4, fid);

fread(&ushSamplesInScanTmp, 1, 2, fid);

fread(&ushUsedChannelsTmp, 1, 2, fid);

unsigned short ushLine;
fread(&ushLine, 1, 2, fid);

unsigned short ushAcquisition;
fread(&ushAcquisition, 1, 2, fid);

unsigned short ushSlice;
fread(&ushSlice, 1, 2, fid);

unsigned short ushPartition;
fread(&ushPartition, 1, 2, fid);

unsigned short ushEcho;
fread(&ushEcho, 1, 2, fid);

unsigned short ushPhase;
fread(&ushPhase, 1, 2, fid);
unsigned short ushRepetition;
fread(&ushRepetition, 1, 2, fid);

fseek(fid, 18+8+4+2+2+16+28, SEEK_CUR); 


unsigned short ushChannelID;
fread(&ushChannelID, 1, 2, fid);
fseek(fid, 2, SEEK_CUR);



if (oldMostSig == -1)
{
	ushSamplesInScan = ushSamplesInScanTmp;
	ushUsedChannels = ushUsedChannelsTmp;

	if (fsize < 0)
		fsize = 6000000000L;

	if (fsize > 0)
	{
		if (allchan)
			totalLine = fsize / (128 + ushSamplesInScan * 4 * 2);
		else
			totalLine = fsize / (128 + ushSamplesInScan * 4 * 2) / ushUsedChannels;
	}
	long int msz = 2 * ushSamplesInScan * 4 * totalLine;  // this does not work when size is too big
	
		data = (float*)malloc(msz);

		phInd = (unsigned short *)malloc(2 * totalLine);
		parInd = (unsigned short *)malloc(2 * totalLine);
		Rep = (unsigned short *)malloc(2* totalLine);
		slcInd = (unsigned short *)malloc(2 * totalLine);
		echoInd = (unsigned short *)malloc(2 * totalLine);
		chInd = (unsigned short *)malloc(2 * totalLine);


		if (data == NULL || phInd == NULL || parInd == NULL || Rep == NULL || slcInd == NULL || echoInd == NULL || chInd == NULL)
		{

			fclose(fid);
			printf("Error allocate memory for the data\n");
			return;
		}
}

if (oldMostSig != aulEvalInfoMaskMostSig)
{
	oldMostSig = aulEvalInfoMaskMostSig;

	short str[32];
	dec2bin(aulEvalInfoMaskMostSig, 32, str);
	if (str[0] == 1) acqEnded = true;
	rawdataCor=(str[10]==1)?true:false;
}


if (acqEnded) break;

if (mdhCounter==totalLine)
{
    printf("Break before all lines are read\n");
    break;
}
if (!allchan && ushChannelID != chan)
 fseek(fid, 2 * ushSamplesInScan * 4, SEEK_CUR);
else
{
   fread(&data[mdhCounter*2 * ushSamplesInScan], 2 * ushSamplesInScan,4, fid);
   
   if (rawdataCor&&nrhs>2)
   {
     for (int i=0;i<ushSamplesInScan;i++)
     {

      float rl=data[mdhCounter*2 * ushSamplesInScan+i*2];
      float im=data[mdhCounter*2 * ushSamplesInScan+i*2+1];

        data[mdhCounter*2 * ushSamplesInScan+i*2]=rl*sclr[ushChannelID]-im*scli[ushChannelID];
        data[mdhCounter*2 * ushSamplesInScan+i*2+1]=rl*scli[ushChannelID]+im*sclr[ushChannelID];

       
     }
       
   }
       
   if ((allchan && ushChannelID == 0) || !allchan)
   {
       int cnt2;
       if (allchan)
	     cnt2 = mdhCounter / ushUsedChannels;
       else
         cnt2 = mdhCounter;
       
	   phInd[cnt2] = ushLine;
	   parInd[cnt2] = ushPartition;
	   Rep[cnt2] = ushRepetition;
	   slcInd[cnt2] = ushSlice;
	   echoInd[cnt2] = ushEcho;
	   chInd[cnt2] = ushChannelID;
   }

	mdhCounter++;
	if (mdhCounter % 10000 == 0)
	{
		printf("% d \n", mdhCounter);
		mexEvalString("drawnow;");
	}
	
}


}

fclose(fid);

if (nlhs>0)
{

    
	int odims[] = { ushSamplesInScan, mdhCounter };
	plhs[0] = mxCreateNumericArray(2, odims, mxSINGLE_CLASS, mxCOMPLEX);
	if (plhs[0] == NULL) 
	{
		printf("cannot create output array\n");
		return;
	}

	float *datar = (float *)malloc(mdhCounter * 4 * ushSamplesInScan);
	float *datai = (float *)malloc(mdhCounter * 4 * ushSamplesInScan);

	for (int ii = 0; ii < mdhCounter * ushSamplesInScan; ii++)
	{
		datar[ii] = data[ii * 2];
		datai[ii] = data[ii * 2 + 1];
	}

	memcpy(mxGetPr(plhs[0]), datar, mdhCounter* sizeof(float)  * ushSamplesInScan);
	memcpy(mxGetPi(plhs[0]), datai, mdhCounter* sizeof(float)  * ushSamplesInScan);

	free(datar);
	free(datai);
     
    /*
    int odims[] = { 2,ushSamplesInScan, mdhCounter };
	plhs[0] = mxCreateNumericArray(3, odims, mxSINGLE_CLASS, mxREAL);
	if (plhs[0] == NULL) 
	{
		printf("cannot create output array\n");
		return;
	}
	memcpy(mxGetPr(plhs[0]), data, 2*mdhCounter* sizeof(float)  * ushSamplesInScan);
	*/
}

int cnt2;
if (allchan)
cnt2 = mdhCounter / ushUsedChannels;
else
cnt2 = mdhCounter;

if (nlhs > 1)
{
	int odims[] = { 1, cnt2 };
	plhs[1] = mxCreateNumericArray(2, odims, mxUINT16_CLASS, mxREAL);
	if (plhs[1] == NULL) /* check that mem was allocated */
	{
		printf("cannot create output array\n");
		return;
	}
	memcpy(mxGetPr(plhs[1]), phInd, cnt2* 2);
}


if (nlhs > 2)
{
	int odims[] = { 1, cnt2 };
	plhs[2] = mxCreateNumericArray(2, odims, mxUINT16_CLASS, mxREAL);
	if (plhs[2] == NULL) /* check that mem was allocated */
	{
		printf("cannot create output array\n");
		return;
	}
	memcpy(mxGetPr(plhs[2]), parInd, cnt2* 2);
}

free(data);

free(phInd);
free(parInd);
free(Rep);
free(slcInd);
free(echoInd);
free(chInd);

}



/*
bits = num2str(dec2bin(evalInfo));
bits = fliplr([repmat('0',1,32-length(bits)),bits]);
setFlags = find(bits=='1')-1;

mdhBitFields.MDH_ACQEND = any(setFlags==0);
mdhBitFields.MDH_RTFEEDBACK = any(setFlags==1);
mdhBitFields.MDH_HPFEEDBACK = any(setFlags==2);
mdhBitFields.MDH_ONLINE    = any(setFlags==3);
mdhBitFields.MDH_OFFLINE   = any(setFlags==4);
mdhBitFields.MDH_LASTSCANINCONCAT = any(setFlags==8);       % Flag for last scan in concatination
mdhBitFields.MDH_RAWDATACORRECTION = any(setFlags==10);      % Correct the rawadata with the rawdata correction factor
mdhBitFields.MDH_LASTSCANINMEAS = any(setFlags==11);      % Flag for last scan in measurement
mdhBitFields.MDH_SCANSCALEFACTOR = any(setFlags==12);      % Flag for scan specific additional scale factor
mdhBitFields.MDH_2NDHADAMARPULSE = any(setFlags==13);      % 2nd RF exitation of HADAMAR
mdhBitFields.MDH_REFPHASESTABSCAN = any(setFlags==14);      % reference phase stabilization scan
mdhBitFields.MDH_PHASESTABSCAN = any(setFlags==15);      % phase stabilization scan
mdhBitFields.MDH_D3FFT     = any(setFlags==16);      % execute 3D FFT
mdhBitFields.MDH_SIGNREV   = any(setFlags==17);      % sign reversal
mdhBitFields.MDH_PHASEFFT  = any(setFlags==18);      % execute phase fft
mdhBitFields.MDH_SWAPPED   = any(setFlags==19);      % swapped phase/readout direction
mdhBitFields.MDH_POSTSHAREDLINE = any(setFlags==20);      % shared line
mdhBitFields.MDH_PHASCOR   = any(setFlags==21);      % phase correction data
mdhBitFields.MDH_PATREFSCAN = any(setFlags==22);      % additonal scan for PAT reference line/partition
mdhBitFields.MDH_PATREFANDIMASCAN = any(setFlags==23);      % additonal scan for PAT reference line/partition that is also used as image scan
mdhBitFields.MDH_REFLECT   = any(setFlags==24);      % reflect line
mdhBitFields.MDH_NOISEADJSCAN = any(setFlags==25);      % noise adjust scan --> Not used in NUM4
mdhBitFields.MDH_SHARENOW  = any(setFlags==26);      % all lines are acquired from the actual and previous e.g. phases
mdhBitFields.MDH_LASTMEASUREDLINE = any(setFlags==27);      % indicates that the current line is the last measured line of all succeeding e.g. phases
mdhBitFields.MDH_FIRSTSCANINSLICE = any(setFlags==28);      % indicates first scan in slice = any(setFlags==needed for time stamps)
mdhBitFields.MDH_LASTSCANINSLICE = any(setFlags==29);      % indicates  last scan in slice = any(setFlags==needed for time stamps)
mdhBitFields.MDH_TREFFECTIVEBEGIN = any(setFlags==30);      % indicates the begin time stamp for TReff = any(setFlags==triggered measurement)
mdhBitFields.MDH_TREFFECTIVEEND = any(setFlags==31);


*/
