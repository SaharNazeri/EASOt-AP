function dsac=readsacs(file_name,inlect,type)
% dsac=readsac(file_name,<inlect>)
% Reading  a seismic record with the SAC format
% dsac=readsac(file_name) with file_name : name of the SAC file (the path of the directory
% containing the file is possible).
% inlect (optional) : type of reading (to open the file in the desired mode (r: reading only,
% r+ reading and writing ... etc)), for more information see the help of the function fopen.
%
% Output :
% dsac is a structure
%
% Modified on March 09, 2004 by A. KAVIANI, to be able to use all the 10 fields
% associated with the picked time of some phases: h1(11:20) and h3(49:128)
% related to the fields T0-T9 and KT0-KT9 of the SAC header.
% See the document of the SAC for more information.

% par defaut on lit les signaux en mode 'r'

% By default signals are read in 'r' mode
if nargin ==1
    typ='pc';
    modlect='r';
elseif nargin ==2
    typ='pc';
    modlect=inlect;
else
    typ=type;
    modlect=inlect;
end
% lecture des signaux
% -------------------

% Suivant le type des signaux la lecture s'effectue de mani�re differente
% pour les stations Minititan :
%	 'l' - IEEE floating point with little-endian byte ordering
% pour les stations Ceispace :
%	'b' - IEEE floating point with big-endian byte ordering

% For the Minititan stations, the reading is done following the signals
% type:
%	 'l' - IEEE floating point with little-endian byte ordering
% for the Ceispace stations :
%	'b' - IEEE floating point with big-endian byte ordering

boolceis=0;	% 0 = Minititan		1 = CEIS
%if strcmp(typ,'sun')
%	fidev=fopen(file_name,modlect,'b');
%end
%if strcmp(typ,'pc')
%	fidev=fopen(file_name,modlect,'l');
%end

fidev=fopen(file_name,modlect,'l');

if fidev > 0

    h1=fread(fidev,70,'float');		% -----------------
    h2=fread(fidev,40,'long');		% reading the header
    h3=fread(fidev,192,'uchar');	% -----------------

    % to test to know in which order are the bytes
    % (if we are reading a Ceispace or a Minititan)

    if h2(7)~=6
        boolceis=1;
        fclose(fidev);
        fidev=fopen(file_name,modlect,'l');
        h1=fread(fidev,70,'float');
        h2=fread(fidev,40,'long');
        h3=fread(fidev,192,'uchar');
    end
%    h2
%    return
    % reading the trace

    tamp=fread(fidev,inf,'float');

    dsac.h1=h1;
    dsac.h2=h2;
    dsac.h3=h3;

    % reading the parameters of type float
    % ------------------------------------

    dsac.tau=h1(1);		% Sampling interval
    dsac.DEPMIN=h1(2);
    dsac.DEPMAX=h1(3);
    dsac.DEPMEN=h1(5);
    dsac.bt=h1(6);		% Begin time (1st sample shifting)
    dsac.et=h1(7);		% End time
    dsac.evt0=h1(8);	% origin time of the event / reference

    
    dsac.t0=h1(11);
    dsac.t1=h1(12);
    dsac.t2=h1(13);	% Picked phases time (10 possibles)
    dsac.t3=h1(14);
    dsac.t4=h1(15);
    dsac.t5=h1(16);
    dsac.t6=h1(17);
    dsac.t7=h1(18);
    dsac.t8=h1(19);
    dsac.t9=h1(20);


    dsac.slat=h1(32);
    dsac.slon=h1(33);	% station coordinates
    dsac.salt=h1(34);

    dsac.elat=h1(36);
    dsac.elon=h1(37);	% Event coordinates
    dsac.edep=h1(39);

    dsac.par0=h1(41);
    dsac.par1=h1(42);	% Free space (10 possibles)
    dsac.par2=h1(43);
    dsac.par3=h1(44);
    dsac.par4=h1(45);
    dsac.par5=h1(46);
    dsac.par6=h1(47);
    dsac.par7=h1(48);
    dsac.par8=h1(49);
    dsac.par9=h1(50);

    dsac.distkm=h1(51);	% distance ev-sta in km
    dsac.az=h1(52);		% azimuth ev-sta
    dsac.baz=h1(53);	% backazimuth ev-sta
    dsac.dist=h1(54);	% distance ev-sta in degrees
    dsac.cmpaz=h1(58);	% azimuth
    dsac.cmpinc=h1(59);	% incidence angle


    % reading the parametres of type long
    % -----------------------------------

    % First sample origin time
    dsac.an=h2(1);
    dsac.jr=h2(2);
    dsac.hr=h2(3);
    dsac.mn=h2(4);
    dsac.sec=h2(5);
    dsac.msec=h2(6);

    dsac.npts=h2(10);	% Number of samples

    %h2(28)
    dsac.us1=h2(28);%disp(h2(28));
    dsac.us2=h2(29);
    dsac.us3=h2(30);
    dsac.us4=h2(31);	%Free Space
    dsac.us5=h2(32);
    dsac.us6=h2(33);
    dsac.us7=h2(34);
    dsac.us8=h2(35);

    % reading the parameters of type uchar
    % ------------------------------------

    %dsac.staname=char(h3(1:8))';	% station name
    dsac.staname=rm_blanc(h3(1:8))';
    dsac.NET=rm_blanc(h3(169:170))';
    dsac.evname=rm_blanc(h3(9:24))';	% Event name (region)

    dsac.kt0=rm_blanc(h3(49:56))';
    dsac.kt1=rm_blanc(h3(57:64))';
    dsac.kt2=rm_blanc(h3(65:72))';	% Identifying the picked phases h1(11:20)
    dsac.kt3=rm_blanc(h3(73:80))';
    dsac.kt4=rm_blanc(h3(81:88))';
    dsac.kt5=rm_blanc(h3(89:96))';
    dsac.kt6=rm_blanc(h3(97:104))';
    dsac.kt7=rm_blanc(h3(105:112))';
    dsac.kt8=rm_blanc(h3(113:120))';
    dsac.kt9=rm_blanc(h3(121:128))';

    dsac.kus1=rm_blanc(h3(137:144))';
    dsac.kus2=rm_blanc(h3(145:152))';	% Free Space
    dsac.kus3=rm_blanc(h3(153:160))';

    dsac.kcomp=rm_blanc(h3(161:168))';	% Component
    dsac.kinst=rm_blanc(h3(185:192))';	% Instrument name


    % Reading the trace
    % -------------------

    dsac.trace=tamp;

    if strcmp(modlect,'r+')
        dsac.fidev=fidev;
    else
        fclose(fidev);
    end

end
if fidev==-1
    disp('Error in reading - File not found')
    dsac.tau=-1;

end

function ch1=rm_blanc(ch)

% Looking for string blancs  and removing them
if ischar(ch)
    ch1=ch(find(double(ch)~=32));
else
    ch1=ch(find(ch~=32));
    ch1=char(ch1);
end
