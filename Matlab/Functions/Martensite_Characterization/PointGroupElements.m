function [MatElems,PGN]=PointGroupElements(varargin)
% PointGroupElements returns the matrices corresponding to the symmetry
%                    operations in the user-provided 3D point group.
%
% PointGroupName can be provided in either Schoenflies or International
% notation. Matrix elements are returned as a cell structure.
%
% VARIABLE ARGUMENTS and EXAMPLE USAGE:
%  Point group elements from point group name:
%    MatElems=PointGroupElements('m-3m')%
%  ... from Space Group Number:
%    MatElems=PointGroupElements(225,'SpaceGroupNumber')
%  ... from TSL-EDAX Laue Group code:
%    MatElems=PointGroupElements(43,'TSL')
%  ... from Point group number code:
%    MatElems=PointGroupElements(32,'PointGroupNumber')
%
% KNOWN ISSUES:
%  The returned matrices may not be correct for 16 out of the 230 space
%  groups. The number of matrices returned is correct (see OJ Curnow, 
%  J Chem Ed 84 (2007) p1430 for a clever mnemonic), but the generator
%  matrices vary for these in the CTEMsoft code (from which I borrowed
%  heavily, and which is licensed under GPLv2(C) 2001, 2002 Marc De Graef.
%  I found that these produce some differences in the final set of
%  matrices. The important matrices for my work so far are all covered
%  here, so I'm not putting too much time into figuring this out yet.
%  However... my hypothesis is that: 9 of these depend on your choice of
%  axes in the hexagonal and trigonal cases; two of these may be errors
%  (see comments in code); and the remaining 5 may depend on your choice of
%  axes in the tetragonal system. Note that there is no distinction in the
%  generator symmetry elements in reference [3] p226, and it is implied
%  elsewhere (e.g. Cayron, Acta Cryst A62 (2006) p21, sect 3.2)that all
%  point group symmetry elements should be the same. Thus, it remains
%  possible that some error(s) remain in subfunction MakeGenStr.
%
% REQUIRES: TVec2RMat
%
% REFERENCES:
%   [1] M De Graef, Introduction to Conventional Transmission Electron
%       Microscopy. Cambridge (2003)
%   [2] M De Graef, CTEMsoft (Fortran 90 code supplementing [1])
%   [3] G Burns, AM Glazer, Space Groups for Solid State Scientists, 2nd
%       Edition. Boston: Academic Press, 1990.
%   [4] GD Nigam, "A Matrix Approach to Point Group Symmetries," 
%       J Chem Ed 60 (1983) p919.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       THIS FUNCTION IS PART OF THE ORIENTATION VISUALIZER LIBRARY       %
%                                                                         %
%           (C) Copyright 2012 Eric Payton, eric.payton@rub.de            %
%                                                                         %
%  Licensed under the EUPL, Version 1.1 only (the "License"); You may not %
%  use this work except in compliance with the License. You may obtain a  %
%  copy of the License (in 19 languages) at:                              %
%                                                                         %
%                     http://ec.europa.eu/idabc/eupl5                     %
%                                                                         %
%  Unless required by applicable law or agreed to in writing, software    %
%  distributed under the Licence is distributed on an "AS IS" basis,      %
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or        %
%  implied. See the License for the specific language governing           %
%  permissions and limitations under the License.                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Interpret varargins
if length(varargin)==2
    if strcmpi(varargin{2},'TSL')
        PGN=TSL_Laue_Group_Name(varargin{1});
    elseif strcmpi(varargin{2},'SpaceGroupNumber')
        PGN=PG_Name_From_SG_Number(varargin{1});
    elseif strcmpi(varargin{2},'PointGroupNumber')
        PGN=varargin{1};
    end
else
        PGN=InterpretPointGroupName(varargin{1});
end    

%% Interpret the input and get the generator matrices

GenStr=MakeGenStr(PGN);
M=InterpretGenStr(GenStr);
assignin('base','M',M)

% Multiply the generator matrices together until no new matrices are
% created. The correct number of matrices are produced (see OJ Curnow,
% Chemical Education Today 84 (2007) p1430 for a nice mnemonic trick.)
% The exact matrices produced in the present code may, however, not all be
% correct. See "Known Issues" in the header.
n1=1;n2=length(M(:,1))+1;k=n2;
while n1<n2
    for i=n1:n2
        for j=1:n2
            tmp=[M(i,1:3);M(i,4:6);M(i,7:9)]*[M(j,1:3);M(j,4:6);M(j,7:9)];
            M(k,:)=[tmp(1,1:3) tmp(2,1:3) tmp(3,1:3)];
            k=k+1;
        end
    end
    M=unique(M,'rows');
    n1=n2;
    n2=length(M(:,1));
    k=n2+1;
end

% Turn the result into a cell structure of matrices
MatElems=TVec2RMat(M);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PGN=TSL_Laue_Group_Name(TSL_Laue_Group_Name)
Opt={'1','2','20','22','4','42','3','32','6','62','23','43'};
Key=[2 5 5 8 11 15 17 20 23 27 29 32];
x=zeros(length(Key));
for i=1:length(Opt),x(i)=strcmpi(Opt{i},TSL_Laue_Group_Name);end
loc=find(x==1);
if length(loc)~=1,
    error('Invalid Point Group Name.');
else
    PGN=Key(loc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PGN=PG_Name_From_SG_Number(SpaceGroupNumber)
SGN=str2double(SpaceGroupNumber);
switch SGN
    case 1,                 PGN= 1;
    case 2,                 PGN= 2;
    case num2cell(  3:  5), PGN= 3;
    case num2cell(  6:  9), PGN= 4;
    case num2cell( 10: 15), PGN= 5;
    case num2cell( 16: 24), PGN= 6;
    case num2cell( 25: 46), PGN= 7;
    case num2cell( 47: 74), PGN= 8;
    case num2cell( 75: 80), PGN= 9;
    case num2cell( 81: 82), PGN=10;
    case num2cell( 83: 88), PGN=11;
    case num2cell( 89: 98), PGN=12;
    case num2cell( 99:110), PGN=13;
    case num2cell(111:122), PGN=14;
    case num2cell(123:142), PGN=15;        
    case num2cell(143:146), PGN=16;
    case num2cell(147:148), PGN=17;        
    case num2cell(149:155), PGN=18;        
    case num2cell(156:161), PGN=19;
    case num2cell(162:167), PGN=20;
    case num2cell(168:173), PGN=21;
    case 174,               PGN=22;
    case num2cell(175:176), PGN=23;        
    case num2cell(177:182), PGN=24;
    case num2cell(183:186), PGN=25;
    case num2cell(187:190), PGN=26;
    case num2cell(191:194), PGN=27;
    case num2cell(195:199), PGN=28; 
    case num2cell(200:206), PGN=29; 
    case num2cell(207:214), PGN=30;
    case num2cell(215:220), PGN=31; 
    case num2cell(221:230), PGN=32;
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PGN=InterpretPointGroupName(PointGroupName)
% Options include some likely misspellings and both Schoenflies and
% International notations.
Opt={'1' 'C1' ... %1
    '-1' 'S2' 'Ci' ... %2
    '2' 'C2' ... %3
    'm' 'C1h' 'Cs' ... %4
    '2/m' 'C2h' ... %5
    '222' 'D2' 'V'... %6
    'mm2' 'C2v' ... %7
    'mmm' 'D2h' 'Vh' ... %8
    '4' 'C4' ... %9
    '-4' 'S4' ... %10
    '4/m' 'C4h' ... %11
    '422' 'D4' ... %12
    '4mm' 'C4v' ... %13
    '42m' 'D2d' 'Vd' '-42m' ... %14
    '4/mmm' 'D4h' ... %15
    '3' 'C3' ... %16
    '-3' 'S6' 'C31' 'C3i' ... %17
    '32' 'D3' ... %18
    '3m' 'C3v' ... %19
    '-3m' 'D3d' ... %20
    '6' 'C6' ... %21
    '-6' 'C3h'... %22
    '6/m' 'C6h' ... %23
    '622' 'D6' ... %24
    '6mm' 'C6v' ... %25
    '6m2' '-6m2' '-62m' '62m' 'D3h'... %26
    '6/mmm' 'D6h' ... %27
    '23' 'T' ... %28
    'm3' 'm-3' 'Th' ... %29
    '432' 'O' ... %30
    '43m' '-43m' 'Td' ... %31
    'm3m' 'm-3m' 'Oh'}; %32
Key=[1 1 2 2 2 3 3 4 4 4 5 5 6 6 6 7 7 8 8 8 9 9 10 10 11 11 12 12 13 ...
     13 14 14 14 14 15 15 16 16 17 17 17 17 18 18 19 19 20 20 21 21 ...
     22 22 23 23 24 24 25 25 26 26 26 26 26 27 27 28 28 29 29 30 30 ...
     31 31 31 31 32 32 32];
x=zeros(length(Key));
for i=1:length(Opt),x(i)=strcmpi(Opt{i},PointGroupName);end
loc=find(x==1);
if length(loc)~=1,
    error('Invalid Point Group Name.');
else
    PGN=Key(loc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M=InterpretGenStr(GenStr)
% Loop through codes in GenStr and get the corresponding generator matrices
n=length(GenStr);
M=zeros(n,9);
for i=1:n
    tmp=MakeGenerator(GenStr(i));
    M(i,:)=[tmp(1,1:3) tmp(2,1:3) tmp(3,1:3)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GenStr=MakeGenStr(PGN) 
% Assign the generator string based on the point group number

%       Code      Point Group     Space Groups     Bravais Lattice    Name  Notes
Opt={   'a'    ...  %  1               1       --- Triclinic           1
        'h'     ... %  2               2                              -1
        'c'     ... %  3              3-5      --- Monoclinic          2
        'j'     ... %  4              6-9                              m
        'ch'    ... %  5             10-15                            2/m
        'bc'    ... %  6             16-24     --- Orthorhombic       222
        'bj'    ... %  7             25-46                            mm2   39: bc in MDG code (why not in 222?)
        'bch'   ... %  8             47-74                            mmm
        'bg'    ... %  9             75-80     --- Tetragonal          4
        'bm'    ... % 10             81-82                            -4
        'bgh'   ... % 11             83-88                            4/m
        'bgc'   ... % 12             89-98                            422   90: bg and no c in MDG code (why not in Tetrag 4?)
        'bgj'   ... % 13             99-110                           4mm
        'bmc'   ... % 14            111-122                          -42m   115-120: bmj
        'bgch'  ... % 15            123-142                          4/mmm
        'n'     ... % 16            143-146    --- Trigonal            3
        'nh'    ... % 17            147-148        (Rhombohedral)     -3
        'ne'    ... % 18            149-155                           32    149,151,153: nf
        'nk'    ... % 19            156-161                           3m    157,159: nl
        'neh'   ... % 20            162-167                          -3m    162,163: nf
        'nb'    ... % 21            168-173    --- Hexagonal          6        
        'ni'    ... % 22              174                            -6
        'nbh'   ... % 23            175-176                          6/m
        'nbe'   ... % 24            177-182                          622
        'nbk'   ... % 25            183-186                          6mm
        'nie'   ... % 26            187-190                         -6m2    187,188: nik
        'nbeh'  ... % 27            191-194                         6/mmm
        'bcd'   ... % 28            195-199    --- Cubic              23
        'bcdh'  ... % 29            200-206                           m3
        'bcde'  ... % 30            207-214                          432
        'bcdl'  ... % 31            215-220                         -43m
        'bcdeh' ... % 32            221-230                         m-3m
        };
    % pass the point group number through the option list
    % to obtain the generator string
    GenStr=Opt{PGN};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gmx=MakeGenerator(t)
% Create the generator matrices. Follows De Graef's approach.
gmx = zeros([3 3]); % preallocate the generator matrix
switch uint8(t)
    case ('a'); gmx(1,1)= 1.0; gmx(2,2)= 1.0; gmx(3,3)= 1.0;
    case ('b'); gmx(1,1)=-1.0; gmx(2,2)=-1.0; gmx(3,3)= 1.0;
    case ('c'); gmx(1,1)=-1.0; gmx(2,2)= 1.0; gmx(3,3)=-1.0;
    case ('d'); gmx(1,3)= 1.0; gmx(2,1)= 1.0; gmx(3,2)= 1.0;
    case ('e'); gmx(1,2)= 1.0; gmx(2,1)= 1.0; gmx(3,3)=-1.0;
    case ('f'); gmx(1,2)=-1.0; gmx(2,1)=-1.0; gmx(3,3)=-1.0;
    case ('g'); gmx(1,2)=-1.0; gmx(2,1)= 1.0; gmx(3,3)= 1.0;
    case ('h'); gmx(1,1)=-1.0; gmx(2,2)=-1.0; gmx(3,3)=-1.0;
    case ('i'); gmx(1,1)= 1.0; gmx(2,2)= 1.0; gmx(3,3)=-1.0;
    case ('j'); gmx(1,1)= 1.0; gmx(2,2)=-1.0; gmx(3,3)= 1.0;
    case ('k'); gmx(1,2)=-1.0; gmx(2,1)=-1.0; gmx(3,3)= 1.0;
    case ('l'); gmx(1,2)= 1.0; gmx(2,1)= 1.0; gmx(3,3)= 1.0;
    case ('m'); gmx(1,2)= 1.0; gmx(2,1)=-1.0; gmx(3,3)=-1.0;
    case ('n'); gmx(1,2)=-1.0; gmx(2,1)= 1.0; gmx(2,2)=-1.0; gmx(3,3)= 1.0;
end % switch case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%