function litzoptEdit()

%  Returns transformer or inductor stranding designs of minimal cost and loss
%
%   LITZOPT  returns a table of AWG wire sizes and
%   strand number and plots these designs relative
%   to the optimal AWG 44 design.  LITZOPT  is for
%   use only with piecewise-linear current waveforms.
%   The user runs LITZOPT to create a data file without
%   transformer chracteristic data.  The user inputs
%   the transformer data %  into that file includingr
%   the output (magnetic field value) of a finite
%   element analysis simulation (e.g. Maxwell).  The
%   user then reruns LITZOPT to optimize
%   the selection of diameter and number of  litz
%   strands considering cost and loss.

%   Version 1.7b updated by R. Scott Mongrain 11SEP19
%   Slight amendments to allow external call and return to main transformer
%   design script.  Additional minor work planned to adjust output.

%   Version 1.7 updated by Jennifer Pollock
%   Updated functions testnan 
%   Updated buildbletable to handle 10 windings

%   Version 1.6.5 updated by Jennifer Pollock
%   Removed bugs so litzopt runs on matlab version 7.2
%   Improved commenting to separate functions and explain some
%   functionality in foreachecitation 

% Version 1.6.4
% Jennifer Pollock  
% Fixed the "max" index bugs that result when the code is run
% in Matlab Version v.7


%Version 1.6.3
%Jennifer Pollock
%   7/8/04 fixed testnan, testnandia, and testnanloss functions to scan the entire matrix for nan

%Version 1.6.2
%   9/3/03 fixed bugs in data entry no legs gapped option.
%   9/3/03 organized data entry in file created to 2d option
%Version 1.6.1

%   8/11/03 added windinght1d function to reduce 1d data entry
%   8/11/03 fixed fea function to generate a new data file - works now
%Version 1.6
%   8/2/03 added the ability plot an orignal design on the optimal design frontier
%   and print cost and loss for original design to 1d, 2d, fea analysis
%   8/2/03 changed 1d data file to take in the bobbin area, 'areabb', set breadth of bobbin
%   window equal to breadth of core window, set bobbin height, h = areabb/bw and set 'a' equal to
%   the bobbin height divided by the number of windings - a = h/numwind
%   7/30/03 added function roundtoAWGrange that takes any AWG vector 
%   and rounds the dia matrix to values in AWG vector
%   7/17/03 jdp switched gap = b or c in B2avgcalc1d
%
%   Version 1.5
%   Jennifer Pollock
%   updated the local version with slope matching optimization 6/22/03
%
%   Version 1.4.1 5/1/03
%   Charlie Sullivan
%   Corrected comments on 1.4 and added semicolons where missing
%   Corrected gap interpretation to do all legs gapped for 'A...'  Removed 2D outer legs only, since that isn't supported for 2D. 
%   
%   Version 1.4  4/24/03
%   Charlie Sullivan
%   Note: gaplen is length of each gap, not total.
%   Updates: looked for best of 1.3, 1.29.  Seems to be just cosmetic changes.
%   Added space between matrix elements to avoid annoying error messages
%   line 1396 in B2avgcalc2d, made numwind plus gap = numwind +2 where it was 1, and other changes in ungap
%   Made centergap at left of window.
%   removed length(yzero...) and replaced with correct dimension, named n1yzp.  Line 1540 and the half dozen after that.
%   fixed gap excitation magnitude in foreachexcitation
%   note:  does rectfldcalc for unexcited windings, with zero as amplitude--waste of time for larger number of windings--fixed--see below
%   fixed B^2 = Bx^2 + by^2 not (bx +by)^2
%   fixed Hx, Hy calc to accumulate over different windings, not throw out old data.
%   Fixed imageclac to trim zeros out of list "excwind" before sending to Rectfldcalc
%   Fixed rectfld calc to use an modified excwind that is the same size as the aexc, bexc, etc. variables--one element for each excited winding
%   or image of an excited winding, instead of one element per orginal winding.  This corrects the calc and eliminates an internal loop, thus
%   speeding it up.   Also and did some other little stuff to speed up, like taking 2*pi multiplication out of inner looop.
%   Tested on 30:30 design from original cost/loss paper
%   Note:  imaginary answers typically result from diagonal elements too large from problem wiht field calc.


%   1.29 updates:  at least some dialog box text improved
%   1.2 Updates;  Simualtions now called internal instead of analytic, 
%   Typos for winex to windex in 1d code fixed
%   Units specified in dialogue boxes
%   All legs, center leg, and outer legs options used to describe 
%   gap positions in dialogue boxes
%   Both outer legs case corrected (excitation vector had too many elements)



%     By  Dr. Charles R. Sullivan and Tarek Abdallah and Jennifer Pollock
%                  Dartmouth College
%           www.thayer.darmouth.edu/inductor
% Please forward questions and comments to charles.sullivan@dartmouth.edu


 % The winding loss estimation algorithm employed in
 % this program is introduced in the paper "Winding
 % Loss Calculation with Multiple Windings, Arbitrary
 % Waveforms and Two-Dimensional Field Geometry"
 % by Dr. Charles R. Sullivan from  % IEEE Industry
 % Applications Society Annual Meeting, Oct. 1999, Phoenix,
 % pp. 2093-2099.  (This paper is superceded by the 2001
 % Transactions on Power Electronics article, but the
 % original is still available if you really want
 % it. )
 %
 %% The stranding cost estimation algorithm is from
%   "Cost-Constrained Selection of Strand Wire and
%   Number in a Litz-Wire Transformer Winding." by
%   Dr. Charles  % R. Sullivan IEEE Transactions
%   on Power Electronics, 16 (2), March 2001,
%   pp. 281-8. from "Proceedings of the 1998 IEEE
%   Industry Applications Society
%   % Annual Meeting, 1998,  pp. 907-911.


%  List Of Functions:
%
%
%  % New:  This file creates a new transformer data
%          file with appropriate matrix and vector dimensions,
%          but no data, to be edited by the user.
%
% 
%   Twodimsimulation: A function that is called when a two-dimensional field simulation is
%                     to be performed. The function prompts the user for data to be
%                     stored in a data file. The function executes a succession of commands
%                     to ask the user for various transformer data in a GUI format. The
%                     function concludes by closing the data file and saving it under the
%                     filename and directory of the user's choice.  The data file is edited
%                     by the user and again called when a two-dimensional field simulation
%                     is to be performed to allow for subsequent strand optimization.
%
%  Onedimsimulation: A function that is called when a one-dimensional field simulation is
%                     to be performed. The function prompts the user for data to be
%                     stored in a data file. The function executes a succession of
%                     commands to ask the user for various transformer data in a GUI
%                     format.  The function concludes by closing the data file and
%                     saving it under the filename and directory of the user's choice.
%                     The data file is edited by the user and again called when a
%                     one-dimensional field simulation is to be performed to allow
%                     for subsequent strand optimization.
%
%   FEA: A function that is called when the user requests field data to be externally
%        computed by finite element analysis and inputted into a Litzopt data file.
%        The function prompts the user for data to be stored in a data file. The function executes
%        a succession of commands to ask the user for various transformer data in a GUI
%        format.  The function concludes by closing the data file and saving it under
%        the filename and directory of the user's choice.  The data
%        file is edited by the user to include the results of a fea-based field
%        simulation to allow for subsequent strand optimization.
%
% 
%    Existing:  This file optimizes the stranding design
%               given a file edited by the user
%
%   averagebsquare:  Generates average value of the squared field
%                     with respect to each of the windings due to the
%                     current in all the windings.  One  input is a three-dimensional
%                     array containing square of the sum of squared
%                     fields averaged over winding 'k' (for terms in page 'k')
%                     due to exited windings 'i' and 'j' (the output of a
%                     finite element analysis or analytic field approximation) on the 
%                     main diagonal and upper tringle.   This input has zero 
%                     terms in lower traingle.  Additional input is 'N' (turns number 
%                     in each winding).  The output is also a three-dimensional 
%                     array of products Bi*Bj at each location (i,j) multiplied by
%                     the prodict of turns number Ni*Nj at location (i,j) on each page.
%
%
%   %Costnloss:  Returns vectors of costs, losses, and stranding
%                for each of the calculated designs where each
%                row designates a different design. Inputs are
%                'Irms' (root mean squared current), 'lt' (turn
%                length), 'N' (number of turns), 'awg' (a vector
%                of strand diameters for each of the optimized
%                designs), 'rhoc' (the resistivity of copper at
%                operating temperature), and 'kl' (a parameter
%                defined by time-rate of change of currents and
%                the modified dynamic loss matrix). The outputs
%                are an array called 'cost' (Cost function of litz
%                wire is assumed to be proportional to
%                '(1+(k/(DC^4))+(kk))*n*lt*N' where 'k' is 1.1e-26,
%               'kk' is 2e-9, 'n' is the number of strands in the
%                design in question 'lt' is the length of a single
%                turn, and 'N' is the number of turns.  where, dc
%                is strand diameter, n is strand number, lt is turn
%                length, and N is turn number),  'actual_loss'
%                (the total winding power lost in each design and is
%                 equal to 'Irms^2*lt*N*rhoc*Fr/n/(pi*dc^2/4)' where
%                 all parameters are as defined above and rhoc is the
%                 resistivity of copper at operating temperature),
%                 and 'n' is defined to be the number of strands
%                 implemented in each design.
%
%    %dCostm:   Approximates the derivative of cost function given
%               strand size and the same cost function is used
%               (see Costnloss)
% 
%     Numstrands: Solves for the number of strands in the optimal
%                 design given strand size and ac to dc resistance
%                 ratio (Fr)
% 
%     Awgtomet:  Converts from American Wire Gauge strand size measure
%                to meters
% 
%     COSTm:     Calculates cost basis proportional to the additional
%                cost per unit mass of a given strand size.  Cost is
%                considered with the function given above (see Costnloss)
%
%
%
%     B2avgcalc1d: Solves for magnitude of magnetic field in one spatial dimension.
%                  Inputs are breadth of winding window, gap locations,uo, and
%                  number of windings
%
%
%
%
%
%     B2avgcalc2d: Generated two-dimensional field solution and output is area-integral
%                  of squared field,  givens are transformer geometry such as, gap
%                  locations and size, winding locations and size, number of requested
%                  images and winding divisions.  Calls function foreachexcitation for
%                  specific excitation
%
%
%
%     foreachexcitation: changes winding excitations in sequence to obtain all field
%                        values described in the SFD method.  inputs are winding positions
%                        and locations in window, gap length and number, breadth and height
%                        of the winding window, gap locations, winding divisions for field
%                        computation, vector describing excitation and number
 %                       of images to be computed
%
%     imagecalc: computes all data needed to use method of images about high perbeablilty
  %              core,  the inputs are dimension and location of each of the windings in
  %              the window, breadth and height of the window, number of divisions for
  %              windings, excited winding vector, and number of images to be computed.
  %              The entire system of imaged excited windings are generated.:
%
%     rectfldcalc: takes in the imaged system of excited windings and computes the
%                  area-average squared magnetic field over the cross-section of a
%                  single transformer winding.  The inputs are the location and
%                  dimensions of the imaged system, the location and dimensions
%                  of the winding being averaged over, breadth and height of the
%                  winding window, and the number of divisions of the winding.
%                  The single area-avergae field value is the output.
%      ungap:  to place Ribbons with same length of gaps and magnitude
%              equal and opposite to winding current; then  ungap the core.
%
%   PWLcurrentplot: takes in 'dt' (vector of time segment lengths for current waveform), 
%                   'I' (matrix of current current values), and 'numwind' to create current 
%                    waveform plots
%
%   Imscurrent: computes mean-square current Ims and root mean squared current Irms
%               from I, dt, numtime, numwind
%
%   find_k_ell:     calculate di/dt for each winding and approximate loss by SFD method
%                   from numwind,numtime,I,Ims,dt,N,len,rhoc,B2avg and returns 'kl'
%
%   find_optimal_designs:   finds the designs for wire sized in vector awg, and also the 
%                           set of full-bobbin designs (all the "2" ones)
%                           numwind, kl, Irms, len, N, rhoc are all pretty obvious.
%                           userfp is the maximum packing factor relative to perfect 
%                           square packing of cylinders that the user specifies; used for finding
%                           full-bobbin designs.  Takes in: awg,numwind,kl,Irms,len,N,rhoc,bb,a,insulbuild,userfp
%                           Returns:[totcost,totloss,lossW,ncorrect,fp,totcost2,totloss2,loss2,n2,totcost3,
%                           totloss3,loss3,n3, awground, n2round,totcost2round, totloss2round, dia
% 
%   find_optimal_designs_curvedata:    Takes in awgfine,numwind,kl,Irms,len,N,rhoc,bb,a,insulbuild,userfp.
%                                       Does the same thing as find_optimal_designs, but with awgfine to create 
%                                       data to plot.  Only returns totcostfine,totlossfine,totcost2fine,totloss2fine
%                                       where totcostfine and totlossfine are the hypothetical optimal design curve
%                                       and totcost2fine and totloss2fine are the buikdabke design curve
%
%   designlabels:   creates designlabels for buildable design table
%
%   plotdesignL:    plots the optimal design frontier
%
%   Xsectionplot:   plots cross-section of winding window to confirm location of windings in the window
%                               
%   findtable:      This function is called by find_optimal_designs and find_optimal_design_curvedata.
%                   Finds the lookup table
%
%   numderiv:   This function calculates the numerical deriv by taking x,y or cost, loss and 
%               returning the derivative and it's location.
%
%   fprime:     This function looks up the derivtive in the look up table givena cost and returns the derivative
%
%   fprimeinverse:  This function looks up the normalized cost in the look up table given a derivative
%
%   dlookup:    This function looks up the diameter in the look up table given a cost
%
%   dlookupinverse: %this function looks up the normalized cost cooresponding to a diameter(awg) in the look up table.
%
%   ff:     This the forward function that given a noralized cost, it looksup the normalized loss
%
%   stranding:  This function takes in Irms,lt,kl,N,rhoc,d2,lossW2 and returns 'n', the number of strands
%               This function calls numstrand
%
%   totaldesignsAWG:    This function is called in find_optimal_designs.  It takes in dia, awg, numwind 
%                       and returns rounded designs 'diaround' matrix of round diameters.
%                        this function will determine all the rounded designs
%                       the criterion are for numbers with a remainer < 0.2, it is rounded down.  
%                       If the remainder falls between 0.2 and 0.8, it is rounded both ways
%                       if the remainder is >0.8 the number id rounded up
%
%   optpackfactor:  This function will determine the packing factor of optimal the optimal 
%                   design considering the insulation build.  It returns the matirx 'fp' and
%                   takes in awg, N, n, bb, a, insulbuild
%
%   buldabledesigns:    This function will determine the a buildable design considering the insulation build
%                       It returns the matrix 'n2' of buildable design stranding
%                       It takes in awg, a, n, fp, N, userfp, bb, insulbuild
%                   
%   totaldesigns:    This function will determine all the rounded designs for rounding 
%                   the stranding.  The criterion are for numbers with a remainer < 0.2, 
%                   it is rounded down.  If the remainder falls between 0.2 and 0.8, it 
%                   is rounded both ways. if the remainder is >0.8 the number is rounded up
%                   It returns awground, n2round and takes in n2, awg, numwind
%
%   singlebuild:    This function deterimes wire diameter with single build insulation
%
%   heavybuild:     This function determines wire diameter with heavy build insulation
%
%   smartround:     This function allows a number or vector to be rounded to a specified 
%                   number of significant digits
%
%   costnloss2:     This is a modified version of Costnloss to determine the cost and 
%                   losses for buildable designs the method remains much the same, 
%                   accept the n2 matrix is passed in to the function and it returns 
%                   just the cost2 and actual_loss2 vectors.
%
%   fullbobbindesigns:  this function will determine the full-bobbin designs using 
%                       the insulation build  
%
%   optdesigntable:     This function generates html table for optimal design given
%                       optdesignlabels, totcost, totloss, numwind, n, dia, fp
%
%   buildabletable:     This function generates Html table for buildable design
%                       uses the smartround function to round cost and loss columns
%                       It takes in: bbdesignlabels, totcost2round, totloss2round, 
%                       numwind,n2round,awground
%
%   rounddesigns:   This function determines how to round designs and number of 
%                   designs to create.  It returns nrounded, awground and takes in awg,n2, numwind   
%
%   testnan:    This function prunes matrix that contains NaN that result from 
%               trying to extrapolate from the lookup table/  It returns inputvector, kcostmatrix, klossmatrix
%               and takes in inputvector, kcostmatrix, klossmatrix
%
%   testnanloss:    %This function prunes matrix that contains NaN that result 
%                   from trying to extrapolate from the lookup table.  It returns:
%                   inputvector, costW, kcostmatrix, klossmatrix and takes in:
%                   inputvector, costW, kcostmatrix, klossmatrix)
%
%   testnandia:     %This function prunes matrix that contains NaN that result from 
%                   trying to extrapolate from the lookup table.  It takes in input vector, 
%                   costW, lossW, kcostmatrix, klossmatrix. It returns: inputvector, costW, 
%                   lossW, kcostmatrix, klossmatrix
%
%   roundtoAWGrange:   this function rounded the matrix of wire diameters to the nearest 
%                       value in awg vector passed in.  It returns the matrix of rounded 
%                       wire sizes.
%
%   originaldesign:     this function determines the cost and loss for the orignal design
%                       if data is entered into the fields 'orig_numstrands' and 'orig_wiresize'.  
%                       It returns the total relative cost and loss for original design
%                       
%   windinght1d:        this function calculates h, ht of core window, hb, ht of bobbin window
%                       bb, breadth of bobbin window and a, ht of each winding using N-number of turns
%                       Irms - RMS currnet, bw - core breadth and areabb-area of bobbin window
%
%
%%%%%%%%%%%%%Start of Main Program%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
answer=questdlg('Select Project File:','Litzopt','New Project: Internal Field Simulation','New Project: External Field Simulation','Existing Project','New Project: Internal Field Simulation');

% Decide if Existing file or New,

if answer(1)=='E'
   existing
elseif answer(14)=='I'
   new
else
    fea
end


%%%%%%%%%%%%%End of Main Program%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new()

%this function setup a new data file for either a 1D or 2D solution

    answer=questdlg('Geometry Considered:','Litzopt','1-D Field Geometery Considered','2-D Field Geometery Considered','2-D Field Geometery Considered');
   
    if answer(1)=='1'
       onedimsimulation
    else
       twodimsimulation
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of new%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

%%%%%%%%%%%%%%%% twodsimulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is run when a new data file is to be created and the 
% user selects a 2 dimensional simulation.   For the analytic sol'n, the input parameters no longer
% include winding volumes and simulated field values as they are calculated
% in another function.  For the 2d case, extra user inputs ars gaplen, xzero, yzero (winding loc's),
% xdiv,ydiv (core divisions)

function twodimsimulation()

prompt={'Number of windings:','Number of time segments for piecewise-linear current waveforms:'};
answer=inputdlg(prompt);
vect=str2double(answer);  
numtime=vect(2);
numwind=vect(1);
uo=4e-7*pi;

answer1=questdlg('Are there any gaps in the core? ','Choose','No legs gapped','Center leg gapped','All legs gapped','No legs gapped');

if answer1(1)=='C'
   gaploc='Center';
   
elseif answer1(1) == 'A'
    gaploc ='All';
else
    gaploc ='No Gap';
    gaplen = 0;
end

answer2=questdlg('What build is the wire insulation?','Choose','Single','Heavy','Single');

if answer2(1)=='S'      %this is single build insulation
   insulbuild='s';
   
else    %this is heavy insulation
    insulbuild = 'h';
end


% check to see if '.m' follows filename

[filename,pathname]=uiputfile('Untitled', 'Save new file as:');
both= and (filename(length(filename)-1)=='.',filename(length(filename))=='m');


% if '.m' does not follow filename, place it there 


if both==0;  % if false statements
   filename=[filename '.m'];
end


% store original pathname and set path to new pathname  

origdirect=cd;
cd(pathname);


% create new file in the desired path

[fid,errormsg]=fopen(filename,'w'); % %#ok<*MSNU>


% begin writing file 

startstr='%%For an example and help, go to www.thayer.dartmouth.edu/inductor/litzopt, or refer to winding data for use with litzopt. \n \n';
startstr=[startstr '%%Please replace the question marks with the pertinent data and leave all else alone....\n \n'];
startstr=[startstr '%%There are ' num2str(numwind) ' windings in your transformer. \n \n'];
startstr=[startstr '%%Resistivity (1.77e-8 is for room temperature copper) \n '];
startstr=[startstr ' rhoc = 1.77e-8; %% [ohm-meters]\n \n'];
startstr=[startstr '%%Perform optimization on designs with wire sizes ranging from American Wire Gauge strand size (32 to 48 is the default) \n'];
startstr=[startstr '  awg = [32:2:48]; %% \n \n'];
startstr=[startstr '%%Number of images of the Winding Window Geometry to be computed in the horizontal direction (5 images is the default) \n'];
startstr=[startstr '  Ximage = 5; %% \n \n'];
startstr=[startstr '%%Number of images of the Winding Window Geometry to be computed in the vertical direction (5 images is the default) \n'];
startstr=[startstr '  Yimage = 5; %% \n \n'];
startstr=[startstr '%%Number of horizontal divisions of a winding (20 divisions is the default) \n'];
startstr=[startstr '  xdiv = 20; %% \n \n'];
startstr=[startstr '%%Number of vertical divisions of a winding (20 divisions is the default) \n'];
startstr=[startstr '  ydiv = 20; %% \n \n'];
startstr=[startstr '%%Breadth of the Core Window in meters \n'];
startstr=[startstr '  bw = ?? ; %% \n \n'];
startstr=[startstr '%%Height of the Core Window in meters \n'];
startstr=[startstr '  h = ?? ; %% \n \n'];
startstr=[startstr '%%Breadth of the Bobbin Window in meters \n'];
startstr=[startstr '  bb = ?? ; %% \n \n'];
startstr=[startstr '%%Height of the Bobbin Window in meters \n'];
startstr=[startstr '  hb = ?? ; %% \n \n'];
startstr=[startstr '%%Length of the gap in meters \n'];
startstr=[startstr '  gaplen = ?? ; %% \n \n'];
startstr=[startstr '%%Maximum Achieveable Packing Factor \n'];
startstr=[startstr '  userfp = ?? ; %% \n \n'];
fprintf(fid,startstr);

%%% Length, create vector of winding lengths 


istr=['%%Average turn length for each winding... ' num2str(numwind) ' elements are required in this vector. \n'] ; 
istr=[istr '  len =[ '];

%%%% loop through windings and designate vector element for each winding length

for i=1:numwind
    istr=[istr ' ?? ']; 
end

% close vector and write it to file 

istr=[istr '];  %%[meters] \n \n'];
fprintf(fid,istr);


%%%  Time Duration    %%%

% create vector of time durations for each segment


dtstr=['%%Duration of each time segment... ' num2str(numtime) ' elements are required in this vector. \n'];
dtstr=[dtstr '  dt =[ '];


% loop through number of time segments and designate a vector element for additional
% time segment duration

for j=1:numtime   
    dtstr=[dtstr ' ?? ']; %#ok<*AGROW>
end

% close vector and write it to file 

dtstr=[dtstr ']; %%[seconds] \n \n'];
fprintf(fid,dtstr);

%%%%% Input Winding location vectors in x and y coords

astr=[ '%%Height of Winding Cross Section for each winding (in h direction)... ', num2str(numwind) ' elements are required in this vector. \n  a=['] ;

% loop through the windings and designate vector element for position in
% each additional winding

for i=1:numwind
    astr=[astr ' ?? '];
end

% close vector and write it to file 

astr=[astr '];  %%[meters] \n \n'];
fprintf(fid,astr);

 bstr=[ '%%%%Width of Winding Cross Section for each winding (in bw direction) ... ', num2str(numwind) ' elements are required in this vector. \n  b=['] ;

% loop through the windings and designate vector element for position in
% each additional winding


for i=1:numwind
    bstr=[bstr ' ?? '];
end

% close vector and write it to file 

bstr=[bstr '];  %%[meters] \n \n'];
fprintf(fid,bstr);

%%%%% Input Winding Starting vectors in x and y coords

xzerostr=[ '%% Location of Center of Winding Cross Section for each winding (in h direction) ... ' num2str(numwind) ' elements are required in this vector. \n  xzero=['] ;

% loop through the windings and designate vector element for position in
% each additional winding

for i=1:numwind
    xzerostr=[xzerostr ' ?? '];
end

% close vector and write it to file 

xzerostr=[xzerostr '];  %%[meters] \n \n'];
fprintf(fid,xzerostr);


yzerostr=[ '%% Location of Center of Winding Cross Section for each winding (in bw direction) ... ' num2str(numwind) ' elements are required in this vector. \n  yzero=['] ;

% loop through the windings and designate vector element for position in
% each additional winding

for i=1:numwind
    yzerostr=[yzerostr ' ?? '];
end

% close vector and write it to file 

yzerostr=[yzerostr '];  %%[meters] \n \n'];
fprintf(fid,yzerostr);

%%%%  I (current)
% loop through the windings and create matrix row to designate each additional winding

for i=1:numwind
    string='%%Current at the end of the last interval is assumed equal to the first current value. \n \n';
    string=[string ' %%Current values for winding ' num2str(numwind-((numwind-i))) ' at the beginning of each time segment... ' num2str(numtime) ' elements are required in this vector. \n   I(' num2str(numwind-((numwind-i))) ',:)= ['];

    % loop through time segments and designate vector element for each additional
    % segment



    for k=1:numtime
        string=[string ' ?? '];
    end

   
   % close vector and write it to file 


   string=[string '];  %%[amperes] \n  \n' ];
   fprintf(fid,string);

end


%%%%  N (turns)

   Nstr=[ '%%Number of turns in each winding... ' num2str(numwind) ' elements are required in this vector. \n  N=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
    Nstr=[Nstr ' ?? '];
end

% close vector and write it to file 

Nstr=[Nstr '];  %%[turns] \n \n'];
fprintf(fid,Nstr);

%%%%  Plot orignal design - create vector of wire diameters

   origdiastr =[ '%%Wire gauge for each winding... ' num2str(numwind) ' elements are required in this vector. \n  orig_wiresize=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
    origdiastr=[origdiastr ' 0 '];
end

% close vector and write it to file 

origdiastr=[origdiastr '];  %%[AWG gauge] \n \n'];
fprintf(fid,origdiastr);

%%%%  Plot orignal design - create vector of number of strands in each winding

   orignumstrandstr =[ '%%Number of strands in each winding... ' num2str(numwind) ' elements are required in this vector. \n  orig_numstrands=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
   orignumstrandstr=[orignumstrandstr ' 0 '];
end

% close vector and write it to file 

orignumstrandstr=[orignumstrandstr '];  %%[number of strands] \n \n'];
fprintf(fid,orignumstrandstr);

% make matlab spit out the gap/core dimensions in the file just as it would display it ordinarily
insulbuildstr=[ '\n %%The wire insulation build is \n  insulbuild= ''' insulbuild '''; \n'];
fprintf(fid,insulbuildstr);   %  print insulation build

gaplocstr=[ '\n %%Location of Gap \n  gaploc= ''' gaploc ''' ; \n'];
% gaplocstr=[ '%%Breadth of Winding Window [meters] \n  bw='  num2str(bw) '; \n \n %%Height of Winding Window \n  h=' num2str(h) '; \n 
%%Breadth of Bobbin Window [meters] \n  bb='  num2str(bb) '; \n \n %%Height of Bobbin Window \n  hb=' num2str(hb) '; \n \n
%%Gaploc \n  gaploc=''' gaploc '''; \n' gaplenstr];
fprintf(fid,gaplocstr);   %  print gap location

%this gives the user tha ablitity to edit datafile immediately
fclose(fid);
msg1=['You have created file  "' (filename) '" in path "' (pathname) '"   Please open this file and enter your data, then rerun litzopt and select "Existing" to calculate the optimal designs for your transformer.' ];
answer=questdlg(msg1,'Choose','Edit Now', 'Ok, I will edit it later', 'default');


if answer(1)=='E'
   edit(filename);
else
   cd(origdirect);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end twodsimulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% onedsimulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is run when a new data file is to be created and the 
% user selects a 1 dimensional simulation.    For the analytic sol'n, the input parameters no longer
% include winding volumes and simulated field values as they are calculated
% in another function


function onedimsimulation()
prompt={'Number of windings:','Number of time segments for piecewise-linear current waveforms:'};
answer=inputdlg(prompt);
vect=str2double(answer);  
numtime=vect(2);
numwind=vect(1);
uo=4e-7*pi; %#ok<*NASGU>

answer1=questdlg('Are there any gaps in the core? ','Choose','No legs gapped','Center leg gapped','All legs gapped','No legs gapped');

if answer1(1)=='C'
   gaploc='Center';
   
elseif answer(1) == 'A'
    gaploc = 'All';
else
    gaploc='No Gap';
end

answer2=questdlg('What build is the wire insulation?','Choose','Single','Heavy','Single');

if answer2(1)=='S'      %this is single build insulation
   insulbuild='s';
   
else    %this is heavy insulation
    insulbuild = 'h';
end

% check to see if '.m' follows filename

[filename,pathname]=uiputfile('Untitled', 'Save new file as:');
both= and (filename(length(filename)-1)=='.',filename(length(filename))=='m');

% if '.m' does not follow filename, place it there 

if both==0;  % if false statements
   filename=[filename '.m'];
end


% store original pathname and set path to new pathname  

origdirect=cd;
cd(pathname);


% create new file in the desired path

[fid,errormsg]=fopen(filename,'w');


% begin writing file 

startstr='%%For an example and help, go to www.thayer.dartmouth.edu/inductor/litzopt, or refer to winding data for use with litzopt. \n \n';
startstr=[startstr '%%Please replace the question marks with the pertinent data and leave all else alone....\n \n'];
startstr=[startstr '%%There are ' num2str(numwind) ' windings in your transformer. \n \n'];
startstr=[startstr '%%Resistivity (1.77e-8 is for room temperature copper) \n '];
startstr=[startstr ' rhoc = 1.77e-8; %% [ohm-meters]\n \n'];
startstr=[startstr '%%Perform optimization on designs with wire sizes ranging from American Wire Gauge strand size (32 to 48 is the default). \n'];
startstr=[startstr '  awg = [32:2:48]; %% \n \n'];
startstr=[startstr '%%Breadth of the Core Window in meters \n'];
startstr=[startstr '  bw = ?? ; %% \n \n'];
startstr=[startstr '%%Area of the Bobbin Window in meters \n'];
startstr=[startstr '  areabb = ?? ; %% \n \n'];
startstr=[startstr '%%Maximum Achieveable Packing Factor \n'];
startstr=[startstr '  userfp = ?? ; %% \n \n'];
fprintf(fid,startstr);

%%% Length, create vector of winding lengths 

istr=['%%Average turn length for each winding... ' num2str(numwind) ' elements are required in this vector. \n'] ; 
istr=[istr '  len =[ '];

%%%% loop through windings and designate vector element for each winding length

for i=1:numwind
    istr=[istr ' ?? '];
end

% close vector and write it to file 

istr=[istr '];  %%[meters] \n \n'];
fprintf(fid,istr);

%%%  Time Duration    %%%

% create vector of time durations for each segment

dtstr=['%%Duration of each time segment... ' num2str(numtime) ' elements are required in this vector. \n'];
dtstr=[dtstr '  dt =[ '];

% loop through number of time segments and designate a vector element for additional
% time segment duration

for j=1:numtime   
    dtstr=[dtstr ' ?? '];
end

% close vector and write it to file 

dtstr=[dtstr ']; %%[seconds] \n \n'];
fprintf(fid,dtstr);

%%%%  I (current)
% loop through the windings and create matrix row to designate each additional winding

for i=1:numwind
    string='%%Current at the end of the last interval is assumed equal to the first current value. \n \n';
    string=[string ' %%Current values for winding ' num2str(numwind-((numwind-i))) ' at the beginning of each time segment... ' num2str(numtime) ' elements are required in this vector. \n   I(' num2str(numwind-((numwind-i))) ',:)= ['];

    % loop through time segments and designate vector element for each additional
    % segment

    for k=1:numtime
        string=[string ' ?? '];
    end

   % close vector and write it to file 

   string=[string '];  %%[amperes] \n  \n' ];
   fprintf(fid,string);

end

%%%%  N (turns)

   Nstr=[ '%%Number of turns in each winding... ' num2str(numwind) ' elements are required in this vector. \n  N=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
    Nstr=[Nstr ' ?? '];
end

% close vector and write it to file 

Nstr=[Nstr '];  %%[turns] \n \n'];
fprintf(fid,Nstr);

%%%%  Plot orignal design - create vector of wire diameters

   origdiastr =[ '%%Wire gauge for each winding... ' num2str(numwind) ' elements are required in this vector. \n  orig_wiresize=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
    origdiastr=[origdiastr ' 0 '];
end

% close vector and write it to file 

origdiastr=[origdiastr '];  %%[AWG gauge] \n \n'];
fprintf(fid,origdiastr);

%%%%  Plot orignal design - create vector of number of strands in each winding

   orignumstrandstr =[ '%%Number of strands in each winding... ' num2str(numwind) ' elements are required in this vector. \n  orig_numstrands=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
   orignumstrandstr=[orignumstrandstr ' 0 '];
end

% close vector and write it to file 

orignumstrandstr=[orignumstrandstr '];  %%[number of strands] \n \n'];
fprintf(fid,orignumstrandstr);


gaplocstr=[ '%%Gap \n  gap=''' gaploc '''; \n'];
fprintf(fid,gaplocstr);   %  print gaplocation

insulationstr=[ '%%The wire insulation build is \n  insulbuild=''' insulbuild '''; \n'];
fprintf(fid,insulationstr);   %  print wire insulation build


fclose(fid);
msg1=['You have created file  "' (filename) '" in path "' (pathname) '"   Please open this file and enter your data, then rerun litzopt and select "Existing" to calculate the optimal designs for your transformer.' ];
answer=questdlg(msg1,'Choose','Edit Now', 'Ok, I will edit it later', 'default');


if answer(1)=='E'
   edit(filename);
else
   cd(origdirect);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%end onedimsimulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% start FEA%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is run whenever a new data file is to be created and the
% user specifies numerical field values are to be used.  The output is a
% data file with the valid matrix and vector dimensions  and blanks 
% (or ??'s) for the user to fill in.

function fea()

prompt={'Number of windings:','Number of time segments for piecewise-linear current waveforms:'};
answer=inputdlg(prompt);
vect=str2double(answer);  
numtime=vect(2);
numwind=vect(1);

answer1=questdlg('Are there any gaps in the core? ','Choose','No legs gapped','Center leg gapped','All legs gapped','No legs gapped');

if answer1(1)=='C'
   gaploc='Center';
   
elseif answer(1) == 'A'
    gaploc = 'All';
else
    gaploc='No Gap';
    gaplen = 0;
end

answer2=questdlg('What build is the wire insulation?','Choose','Single','Heavy','Single');

if answer2(1)=='S'      %this is single build insulation
   insulbuild='s';
   
else    %this is heavy insulation
    insulbuild = 'h';
end

% check to see if '.m' follows filename

[filename,pathname]=uiputfile('Untitled', 'Save new file as:');
both= and (filename(length(filename)-1)=='.',filename(length(filename))=='m');

% if '.m' does not follow filename, place it there 

if both==0;  % if false statements
   filename=[filename '.m'];
end

% store original pathname and set path to new pathname  

origdirect=cd;
cd(pathname);

% create new file in the desired path

[fid,errormsg]=fopen(filename,'w');

% begin writing file 

startstr='%%For an example and help, go to www.thayer.dartmouth.edu/inductor/litzopt, or refer to winding data for use with litzopt. \n \n';
startstr=[startstr '%%Please replace the question marks with the pertinent data and leave all else alone....\n \n'];
startstr=[startstr '%%There are ' num2str(numwind) ' windings in your transformer. \n \n'];
startstr=[startstr '%%Resistivity (1.77e-8 is for room temperature copper) \n '];
startstr=[startstr ' rhoc = 1.77e-8; %% [ohm-meters]\n \n'];
startstr=[startstr '%%Perform optimization on designs with wire sizes ranging from American Wire Gauge strand size (32 to 48 is the default). \n'];
startstr=[startstr '  awg = [32:2:48]; %% \n \n'];
startstr=[startstr '%%Breadth of the Core Window in meters \n'];
startstr=[startstr '  bw = ?? ; %% \n \n'];
startstr=[startstr '%%Height of the Core Window in meters \n'];
startstr=[startstr '  h = ?? ; %% \n \n'];
startstr=[startstr '%%Breadth of the Bobbin Window in meters \n'];
startstr=[startstr '  bb = ?? ; %% \n \n'];
startstr=[startstr '%%Height of the Bobbin Window in meters \n'];
startstr=[startstr '  hb = ?? ; %% \n \n'];
startstr=[startstr '%%Maximum Achieveable Packing Factor \n'];
startstr=[startstr '  userfp = ?? ; %% \n \n'];
fprintf(fid,startstr);



%%% Length, create vector of winding lengths 

istr=['%%Average turn length for each winding...' num2str(numwind) ' elements are required in this vector. \n'] ;
istr=[istr '  len =[ '];


%%%% loop through windings and designate vector element for each winding length

for i=1:numwind

   istr=[istr ' ?? '];

end

% close vector and write it to file 


istr=[istr '];  %%[meters] \n \n'];
fprintf(fid,istr);


%%%  Time Duration    %%%

% create vector of time durations for each segment


dtstr=['%%Duration of each time segment... ' num2str(numtime) ' elements are required in this vector. \n'];
dtstr=[dtstr '  dt =[ '];

% loop through number of time segments and designate a vector element for additional
% time segment duration

for j=1:numtime   
   dtstr=[dtstr ' ?? '];
end

% close vector and write it to file 

dtstr=[dtstr ']; %%[seconds] \n \n'];
fprintf(fid,dtstr);

%%%%  I (current)

% loop through the windings and create matrix row to designate each additional winding

for i=1:numwind

    string='%%Current at the end of the last interval is assumed equal to the first current value. \n \n';
    string=[string ' %%Current values for winding ' num2str(numwind-((numwind-i))) ' at the beginning of each time segment... ' num2str(numtime) ' elements are required in this vector \n  I(' num2str(numwind-((numwind-i))) ',:)= ['];
     
    % loop through time segments and designate vector element for each additional
    % segment

    for k=1:numtime
        string=[string ' ?? '];
    end

   % close vector and write it to file 

   string=[string '];  %%[amperes] \n  \n' ];
   fprintf(fid,string);

end

%%%%  Volume 
  %  create vector of winding volumes

volstr=['%%Volume of each winding... ' num2str(numwind) ' elements are required in this vector. \n  vol=['] ;
  
% loop through the windings and designate vector element for volume of
%  each additional winding

for i=1:numwind
    volstr=[volstr ' ?? '];
end

volstr=[volstr '];  %%[cubic meters] \n \n'];
volstr=[volstr '%%Number of turns in each winding... ' num2str(numwind) ' elements are required in this vector. \n  N=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
    volstr=[volstr ' ?? '];
end

% close vector and write it to file 

volstr=[volstr '];  %%[turns] \n \n'];
fprintf(fid,volstr);

%%%  B (field)

for i= 1:numwind
    bstr=['%%Integral of Bsquared values from simulation with current in winding ' num2str(i) '\n \n'];
    fprintf(fid,bstr);
       
    for j=i:numwind;
       
         if  i~=j
             innerstr= ['%%Integral of B squared values from simulation with current in winding ' num2str(i) ' and ' num2str(j) '\n \n'];
             fprintf(fid,innerstr);
         end
         
         for imat=1:numwind;
             bigstr='       %% Enter integral of B squared ' ;
             bigstr=[bigstr  ' over winding ' num2str(imat) '\n' '       int_B2(' num2str(i)  ',' num2str(j)  ',' num2str(imat) ')= ??;  %%[T^2-m^3 ] \n'];
             fprintf(fid,bigstr);
         end
       
         fprintf(fid,'\n');
     end
 end

 %%%%  Plot orignal design - create vector of wire diameters

   origdiastr =[ '%%Wire gauge for each winding... ' num2str(numwind) ' elements are required in this vector. \n  orig_wiresize=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
    origdiastr=[origdiastr ' 0 '];
end

% close vector and write it to file 

origdiastr=[origdiastr '];  %%[AWG gauge] \n \n'];
fprintf(fid,origdiastr);

%%%%  Plot orignal design - create vector of number of strands in each winding

   orignumstrandstr =[ '%%Number of strands in each winding... ' num2str(numwind) ' elements are required in this vector. \n  orig_numstrands=['] ;

% loop through the windings and designate vector element for turns in
% each additional winding

for i=1:numwind
   orignumstrandstr=[orignumstrandstr ' 0 '];
end

% close vector and write it to file 

orignumstrandstr=[orignumstrandstr '];  %%[number of strands] \n \n'];
fprintf(fid,orignumstrandstr);

%%%%% Input Winding location vectors in x and y coords

astr=[ '%%Height of Winding Cross Section for each winding (in h direction)... ', num2str(numwind) ' elements are required in this vector. \n  a=['] ;

% loop through the windings and designate vector element for position in
% each additional winding

for i=1:numwind
    astr=[astr ' ?? '];
end

% close vector and write it to file 

astr=[astr '];  %%[meters] \n \n'];
fprintf(fid,astr);

bstr=[ '%%%%Width of Winding Cross Section for each winding (in bw direction) ... ', num2str(numwind) ' elements are required in this vector. \n  b=['] ;

% loop through the windings and designate vector element for position in
% each additional winding


for i=1:numwind
    bstr=[bstr ' ?? '];
end

% close vector and write it to file 

bstr=[bstr '];  %%[meters] \n \n'];
fprintf(fid,bstr);

%%%%% Input Winding Starting vectors in x and y coords

xzerostr=[ '%% Location of Center of Winding Cross Section for each winding (in h direction) ... ' num2str(numwind) ' elements are required in this vector. \n  xzero=['] ;

% loop through the windings and designate vector element for position in
% each additional winding

for i=1:numwind
    xzerostr=[xzerostr ' ?? '];
end

% close vector and write it to file 

xzerostr=[xzerostr '];  %%[meters] \n \n'];
fprintf(fid,xzerostr);

yzerostr=[ '%% Location of Center of Winding Cross Section for each winding (in bw direction) ... ' num2str(numwind) ' elements are required in this vector. \n  yzero=['] ;

% loop through the windings and designate vector element for position in
% each additional winding

for i=1:numwind
    yzerostr=[yzerostr ' ?? '];
end

% close vector and write it to file 

yzerostr=[yzerostr '];  %%[meters] \n \n'];
fprintf(fid,yzerostr);

% make matlab spit out the gap location and insulation build in the file just as it would display it ordinarily
insulbuildstr=[ '\n %%The wire insulation build is \n  insulbuild= ''' insulbuild '''; \n'];
fprintf(fid,insulbuildstr);   %  print insulation build

gaplocstr=[ '\n %%Location of Gap \n  gaploc= ''' gaploc ''' ; \n'];
% gaplocstr=[ '%%Breadth of Winding Window [meters] \n  bw='  num2str(bw) '; \n \n 
%%Height of Winding Window \n  h=' num2str(h) '; \n %%Breadth of Bobbin Window [meters] \n  bb='  num2str(bb) '; \n \n 
%%Height of Bobbin Window \n  hb=' num2str(hb) '; \n \n%%Gaploc \n  gaploc=''' gaploc '''; \n' gaplenstr];
fprintf(fid,gaplocstr);   %  print gap location

fclose(fid);

msg1=['You have created file  "' (filename) '" in path "' (pathname) '"   Please open this file and enter your data, then rerun litzopt and select "Existing" to calculate the optimal designs for your transformer.' ];
answer=questdlg(msg1,'Choose','Edit Now', 'Ok, I will edit it later', 'default');

if answer(1)=='E'
   edit(filename);
else
   cd(origdirect);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%end FEA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%  CALCULATIONS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% start of Existing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function imports the data file edited by the user and decides if user-simulated
% field data has been given or the expressions derived from an analytic sol'n are to be used

function existing()

% N by M means N windings and M time segments
% user specifies file to be read

[filename, pathname] = uigetfile('Untitled');

% store original pathname and set path to new pathname

origdirect=cd;
cd(pathname); 

% remove '.m' from filename and execute the file as a string

len=length(filename);
filename=filename(1:1:len-2);
eval(filename);

% Invert matrix dimensions from file so that user can input the vectors in 
% horizontal fashion , N by M means N windings and M time segments

uo=4*pi*1e-7;
N=N'; %#ok<*NODEF>
len=len';
D=size(I);    
orig_numstrands = orig_numstrands';
orig_wiresize = orig_wiresize';
numwind=D(1);   % store the number of windings and time segments
numtime=D(2);
awgfine = min(awg)-4:.1:max(awg)+4;    % Special range for figure and table

%start building new program here

I(:,numtime+1)=I(:,1);    % want currents at beg of interval 1 

%determine the Ims and Irms
[Ims, Irms] = Imscurrent(I, dt, numtime, numwind);

  % if dimfield==2;     % User has selected the 2D sol'n
  % B2avgcalc2d  calculates average squared field data  Inputs are 
  % window breadth gap location, width of gap, uo, number of turns,
  % current, and number of windings
  % Divide by the winding volumes for FEA case prior to runnung hardcalc
           
if sum(strcmp(who,'int_B2')); % User has selected FEA sol'n
   % set verion id to 'f' for fea sol'n
        versionID = 'f';
  % Divide by the winding volumes for FEA case prior to runnung avgbsquare
         for pn=1:numwind
             avg_B2(:,:,pn)=int_B2(:,:,pn)/vol(pn);
         end
 elseif  sum(strcmp(who,'gaplen')); % User has selected 2D sol'n
    % set verion id to 't' for 'two' 2d sol'n
    versionID = 't';
     %%%%%%%%%  B2avgcalc2d  calculates average squared field data 
         [avg_B2]=B2avgcalc2d(a,b,gaplen,gaploc,bw,h,xzero,yzero,xdiv,ydiv,Ximage,Yimage);
    else % % User has selected 1D sol'n
      % set verion id to 'o' for 'one' 1d sol'n
        versionID = 'o';
    
     %%%%%%%%%  B2avgcalc1d  calculates average squared field data 
        [avg_B2]=B2avgcalc1d(bw,gap,uo,numwind);
       
      %calculate bb, hb,h and a from areabb, bw, N, numwind and Irms
        [a,hb,bb, h] = windinght1d(Irms, bw, areabb, N, numwind);
end
     % OK, avg_B2 represents area averages
     % Run averagebsquare in any of the cases for the possible field geometries and get the volume averages 
[B2avg]=averagebsquare(avg_B2,N);
       
% create current waveform figure

PWLcurrentplotL(dt, numwind, I);

%Call a function to do the calculations needed to find k_ell.
kl = find_k_ell(numwind,numtime,I,Ims,dt,N,len,rhoc,B2avg);

%Now call a function to find the actual optimal designs, plus the buildable designs (xxx2 are buildable results), plus...
%special for the paper, the full-bobbin designs (xxx3)
[cost44, totcost,totloss,loss,n,fp,totcost2,totloss2,loss2,n2,totcost3,totloss3,loss3,n3, awground, n2round,totcost2round, totloss2round, dia] = find_optimal_designs(awg,numwind,kl,Irms,len,N,rhoc,bb,a,insulbuild,userfp); %#ok<*ASGLU>

%Do it again for closely spaced wire sizes to get a smooth curve for plotting.
[totcostfine,totlossfine,totcost2fine,totloss2fine] = find_optimal_designs_curvedata(awgfine,numwind,kl,Irms,len,N,rhoc,bb,a,insulbuild,userfp);

% find cost and loss for original component design
[tot_orig_cost, tot_orig_loss] = originaldesign(orig_numstrands, orig_wiresize, Irms, len, N, rhoc, kl, cost44);

%create design labels
optdesignlabels = designlabels(dia);
bbdesignlabels = designlabels(awground);

% %Plot results in design frontiers:

plotdesignsL(tot_orig_cost, tot_orig_loss, bbdesignlabels, optdesignlabels, dia, awground, totloss, totcost, totcost2round, totloss2round, totcost3, totloss3, totcostfine,totlossfine,totcost2fine,totloss2fine)

% %plots optimal design curve with full bobbins curve, hypothetical optimum curve and best buildable curve

%plot winding window cross section
if versionID ~= 'o' 
    XsectionplotL(xzero, yzero, a, b, h, bw, hb, bb);
end

% R. Scott Mongrain - Added a bit of separation and a header to play nice.
fprintf('\nLitzOpt Optimization Report\n')

% displays the optimal design table
htmltable = optdestable(optdesignlabels, totcost, totloss, numwind, n, dia, fp);
disp(htmltable);

% displays the buildable design table
bdtable = buildabletable(bbdesignlabels, totcost2round, totloss2round, numwind,n2round,awground);
disp(bdtable);

% displays original designs cost and loss if an original design was entered
% R. Scott Mongrain - added a bit of separation to play nice.
if tot_orig_cost ~= 0
    title = '\nOriginal Design';
    coststr = ['Relative Cost of Original Design = ' num2str(tot_orig_cost) ''];
    lossstr = ['Loss of Original Design (in Watts) = ' num2str(tot_orig_loss) ''];
    disp(title);
    disp(coststr);
    disp(lossstr);
end

%%%%%%%%%%%%%%%%%%%%%%end of existing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%start of avergebsquare%%%%%%%%%%%%%%%%%%%%%%%%
% Generates spatial average value of the squared field with 
% respect to each of the windings due to the current in all 
% the windings.  Inputs are 'Bsquared',the matrix containing
% square of the sum of area-averaged fields due to windings 'i' and 'j' 
% (the output of a finite element analysis field approximation),
% 'vol' (winding volumes), and 'N' (turns in each winding).  The output
% is a matrix of area-averaged products Bi*Bj at each location (i,j) on each page
% which is the product of fields due to windings 'i' and 'j'.

function [B2avg]=averagebsquare(Bsquared,N)

%determine the number of windings
nw=length(N);
%make a column vector of one with nw rows
one=ones(nw,1);
%make a matrix of the number of turns, i.e. first column is the number of
%turns in primary,Np, second column is the number of turns in the secondary, Ns, etc
Nn=one*N';
%Create a matrix wher Nm(1,1) = Np*Np, Nm(1,2) = Np*Ns, Nm(2,1) = Ns*Np,

Nm=Nn.*Nn';

%replace the (j,i) element with the (i,j) element
for i=1:nw
    for j=1:nw
        Bsquared(j,i,:)=Bsquared(i,j,:);
    end
end

for pn=1:nw
    %take the first page of Bsquared
   B2row=Bsquared(:,:,pn);
   %extracts the diagonal from B2row and makes a row vector 
   Bdiag=spdiags(B2row,0)';
   %multiple Bdiag by one matrix to get a matrix of diagonal elements of
   %B2row
   Bvert=one*Bdiag;
   %obtain the matrix [Bvert(1,1)+Bvert(1,1) Bvert(2,1)+Bvert(1,2); Bvert(1,2)+Bvert(2,1) Bvert(2,2)+Bvert(2,2)]
   allsq.terms=Bvert+Bvert';
   allsq.terms;
   %first, create column vector of diagonal elements in allsq.term, then make that a square matrix with those elements on the diagonal
   %and all other elements zero, then mulitple by 3/2 and a unit matrix of size nw.  finally subtract this matrix from the allsq.terms matrix 
   sq.terms=allsq.terms-eye(nw)*diag(spdiags(allsq.terms,0),0)*(3/2);
   sq.terms;
 
  %take the matrix B2row and subtract sq.terms,, then multiple by the matrix Nm which contains the number of turns squared and divid by 2  
  (B2row-sq.terms); %#ok<*MNEFF>
B2avg(:,:,pn)=(B2row-sq.terms).*Nm/2;
end    
 
%%%%%%%%%%%%%%%%%%%%%%%%%End of avergebsquare%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%start of Costnloss%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns vectors of costs, losses, and stranding for each of the calculated
% designs where each row designates a different design.
% Inputs are 'Irms' (root mean squared
% current), 'lt' (turn length), 'N' (number of turns), 'awg' (a vector of 
% strand diameters for each of the optimized designs), 'rhoc' (the resistivity
% of copper at operating temperature), and 'kl' (a parameter defined by
% time-rate of change of currents and the modified dynamic loss matrix).
% The outputs are an array called 'cost' (Cost function of litz wire is assumed
% to be proportional to '(1+(k/(DC^4))+(kk))*n*lt*N' where 'k' is 1.1e-26, 
% 'kk' is 2e-9, 'n' is the number of strands in the design in question
% 'lt' is the length of a single turn, and 'N' is the number of turns.
% where, dc is strand diameter, n is strand number, lt is turn 
% length, and N is turn number),  'actual_loss' (the total winding power lost
% in each design and is equal to 'Irms^2*lt*N*rhoc*Fr/n/(pi*dc^2/4)' where all 
% parameters are as defined above and rhoc is the resistivity of copper at 
% operating temperature), and 'n' is defined to be the number of strands implemented
% in each design.  


function [cost,actual_loss,n]= costnloss(Irms,lt,kl,N,awg,rhoc)

dc = awgtomet(awg); %convert wire guage to meters for calculation

%we need to calculate the optimal values for Fr at different guages, so do that 
%and put result into Fr[] 


Fr=1.+1./(1.-((2*COSTm(dc))./(dc.*dCOSTm(dc)))); % gives Fr 
%[awg' Fr'];
n = numstrands(dc,kl,Fr);   %numstrands will solve equation (2)

%now n holds a vector of the different number of strands for optimal Fr  at 
%different wire guages 
%now calculate actual_loss 


actual_loss = (Irms^2).*lt.*N.*rhoc.*Fr./n./(pi*dc.^2/4);

% now calculate cost for the cost/loss plot

Co = 0; 
Fr=1.+1./(1.-((2*COSTm(dc))./(dc.*dCOSTm(dc)))); % gives Fr 
cm = COSTm(dc);
cost = (Co+(cm.*(dc.^2).*n)).*(lt*N);

%%%%%%%End of Costnloss%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%startdCostm%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns the derivative of Cm(Dc) (cost).  The input is 'DC' 
% which is the strand diameter

function derivative = dCOSTm(DC)

    k=1.1e-26;
    kk=2e-9;
    derivative=-((6.*k)./(DC.^7))-((2.*kk)./(DC.^3));


%%%%%%%%%%%%%%%%End of dCostm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%start of Numstrands%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Numstrands
% Inputs are 'DC' (strand diameter), 'kl' (a parameter independent of stranding and
% defined by the modified dynamic loss matrix and time rate of change of
% winding currents), and 'Fr' (the ratio of ac to dc resistances).  The output is
% strand number for the design in question

% solves the number of strands for Fr = given values

function numberst = numstrands(DC,kl,Fr)

%Fr=1.+1./(1.-((2*COSTm(DC))./(DC.*dCOSTm(DC)))); % gives Fr 
	Fr1 = Fr-1 ; %our equation uses fr-1 so go ahead and do it
	muo = 4e-7*pi; 
    [numwind,numdesigns]=size(DC);
    klm = kl*ones(1,numdesigns);
	numberst = sqrt(Fr1./klm./(DC.^6)*(4/pi)^3);
    

%%%%%%%%%End of Numstrands%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%start of awgtomet%%%%%%%%%%%%%%%%%%%%%%%
% Converts from American Wire Guage strand size measure to meters. The input is 
% strand diameter as measured by the AWG standard.  The output is the same
% strand diameter in units of meters.

function y = awgtomet(x)

y = 0.0050.*(92).^((36.-x)./39); %convert from awg to inches
y = (2.54e-2)*y;               %convert from inches to meters

%%%%%%%%%%%%%%%%%%%%%%%%%%end of awgtomet%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%  start of  COSTm.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates cost basis proportional to the additional cost per unit mass of a 
% given strand of diameter. The input is 'DC' (the strand diameter) and the 
% output is cost per unit mass assuming a cost function for litz wire
% '1+(1.1e-26/(DC^6))+(2e-9/(DC^2))'. 

function costunitmass = COSTm(DC)

    k = 1.1e-26;
    kk = 2e-9;
    costunitmass = 1+(k./(DC.^6))+(kk./(DC.^2)); % calculate cost for additional litz
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of COSTm %%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of B2avgcalc1d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  B (field) calculation of area-average of B squared over each transformer winding in 1 dimension
%   inputs are breadth of the winding window, gap position, permeability of free space, and number of windings    
%  7/17/03 jdp switched gap = b or c in B2avgcalc1d

function [B2avg]=B2avgcalc1d(bw,gap,uo,numwind)

%the inputs:
    % bw - the breadth of the winding window
    % gap - the gap locations??
    % uo - the permability of free space
    % numwind - the number of windings

%%B=0 is on nongap side
%check to see where gap is if there is one

       %%%%%%  With 1 GAP ON RIGHT %%%%%%%%%%%%%%%%%%%

if gap(1)=='B'

   gapnum=1;   % Run the algorithm for 1 gap on right side
   B1windexcbeg2avg=0;  % spatial integral of squared field before enclosing any current
   B1windexcwind12avg=uo^2/(3*bw^2); % spatial integral of squared field due to excited winding of one winding system
   B1windexcend2avg=uo^2/bw^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   B2windexcbeg2avg=0;
   B2windexcwind12avg=uo^2/(3*bw^2); % spatial integral of squared field due to excited winding of one winding system
   B2windexcmiddle2avg=uo^2/bw^2;
   B2windexcwind22avg=7*uo^2/(3*bw^2);
   B2windsexcend2avg=4*uo^2/(bw^2);  % spatial integral of squared field after enclosing all current
   
elseif gap(1)=='C'
 
    gapnum=2; % Run the algorithm for 1 gap on left side
    B1windexcbeg2avg=uo^2/(bw^2);
    B1windexcwind12avg=uo^2/(3*bw^2); % spatial integral of squared field due to excited winding of one winding system
    B1windexcend2avg=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B2windexcbegavg2=4*uo^2/(bw^2);  % spatial integral of squared field before enclosing any current
    B2windexcwind12avg=7*uo^2/(3*bw^2); % spatial integral of squared field due to excited winding of one winding system
    B2windexcmiddle2avg=uo^2/bw^2;
    B2windexcwind22avg=uo^2/(3*bw^2);
    B2windsexcend2avg=0;  % spatial integral of squared field after enclosing all current
   
else    %Either all or none, and we treat them the same--any excess MMF dropped equally on sides
    
    gapnum=0;  %  Run the algorithm for 0 or 2 gaps
    B1windexcbeg2avg=uo^2/(4*bw^2);  
    B1windexcwind12avg=uo^2/(12*bw^2);
    B1windexcend2avg=uo^2/(4*bw^2);  
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B2windexcbeg2avg=uo^2/bw^2;
    B2windexcwind12avg=uo^2/(3*bw^2);  % spatial integral of squared field due to first excited (current carrying)winding in two winding system
    B2windexcmiddle2avg=0;
    Bmiddle2avg=0;
    B2windexcwind22avg=uo^2/(3*bw^2); %spatial integral of squared field due to second excited (current carrying) winding in two windi
    B2windexcend2avg=uo^2/bw^2;
    
end
          
for pn=1:numwind  % define which winding we average over

  for i= 1:numwind    % USE 1  winding excitation 

        %%%%  Find out where pn is wrt i
        if i<pn % we have a uniform non zero field in our winding,
            B2avg(i,i,pn) =B1windexcend2avg;
        elseif i>pn % we have not enclosed any current and have B(x)=Bbeg in our winding
            B2avg(i,i,pn)=B1windexcbeg2avg;
        else  %i=pn and we have a ramp fcn for B(x)
            B2avg(i,i,pn)=B1windexcwind12avg;
        end
    
  for j=i:numwind;       % try from i+1

      if  i~=j  % USE 2  winding excitation
          
          %%%%%% 2 winding excitation  With GAP %%%%%%%%%%%%%%%%%%%
          %%%  Find out where pn is wrt i and j when j is > then or = to i and i and j are distinct
          allcurrent=and(pn>i,pn>j);
          middle=and(pn>i,pn<j);
          nocurrent=and(pn<i,pn<j);
          
          if allcurrent % we have a uniform non zero field in our winding,
              B2avg(i,j,pn)=B2windexcend2avg;
              B2avg(j,i,pn)=B2avg(i,j,pn);
          elseif middle % we are between excited windings and have B(x)=Bbeg in our winding
              B2avg(i,j,pn)=B2windexcmiddle2avg;
              B2avg(j,i,pn)=B2avg(i,j,pn);
          elseif nocurrent  % we have enclosed no current
              B2avg(i,j,pn)=B2windexcbeg2avg;
              B2avg(j,i,pn)=B2avg(i,j,pn);
          elseif pn==i %i=pn and we have a ramp fcn for B(x)
              B2avg(i,j,pn)=B2windexcwind12avg;
              B2avg(j,i,pn)=B2avg(i,j,pn);
          else  %j=pn and we have a ramp fcn for B(x)
              B2avg(i,j,pn)=B2windexcwind22avg;
              B2avg(j,i,pn)=B2avg(i,j,pn);
          end
%%%%%%   up there %%%  
       end
   end
 end % end i loop
end  % end pn loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of B2avgcalc1d%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of B2avgcalc2d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  B (field) calculation of area-average of B squared over each transformer winding in 2 dimensions
%   inputs are location and dimensions of windings in winding window, gap length, gap location, 
%    breadth and height of the winding window, number of winding divisions to consider the field, and number
%    of images to be computed by method of images about high-permeablity core
% h - height of core window
% bw - breadth of winding window
% a - Height of Winding Cross Section for each winding (in h direction)
% b - Width of Winding Cross Section for each winding (in bw direction)
% xzero - Location of Center of Winding Cross Section for each winding (in h direction)
% yzero - Location of Center of Winding Cross Section for each winding (in bw direction)
% gaplen - length of gaps
% gaploc - the location of the gaps
% xdiv - the number of divisions in the x dimension for field calc
% ydiv - the number of divisions in the y dimension for field calc
% Ximage - the number of images in the x directions for method od images
% Yimage - the number of images in the Y directions for method of images

function [int_B2]=B2avgcalc2d(a,b,gaplen,gaploc,bw,h,xzero,yzero,xdiv,ydiv,Ximage,Yimage)


numwind=length(a);

%check to see where gap is if there is one, compute gap multiplier

if gaploc(1)=='C'  % Center
    gapmult=1;
    numwindplusgap=numwind+1;
    %elseif gaploc(1)=='R'   
    %gapmult=1;
elseif gaploc(1)=='A'  %both center and Outer legs
    gapmult=.5;
    numwindplusgap=numwind+2;
elseif gaploc(1)=='B' %both outer legs, not centergap
    error('Sorry, this gap location option is not supported in this version of the 2D field calculation')
    return %#ok<*UNRCH>
else              %  No gap 
    gapmult=0;
    numwindplusgap=numwind;
    gaplen=0;
end

% Define gapraise

gapraise=bw/2;

% First run 'ungap' to place Ribbons of length gaplen and magnitude 'gapmult' 
% and to ungap the core,  then run 'imagecalc' to obtain field values

       
[agap,bgap,xzerogap,yzerogap]=ungap(gaplen,gapmult,gaploc,gapraise,h);

    
 % the sequence of excitations are called called repeatedely and call 'foreachexcitation' in the pattern of winding excitations
 % begin by initializing excitation vector
 
 excwind=zeros(numwind,1);
 
 
       for i= 1:numwind    % USE 1  winding excitation 
            excwind(i)=1;

            for j=i:numwind;       

                 if  i~=j  % USE 2  winding excitation
                     excwind(i)=1;
                     excwind(j)=1;
                 end
                   
                     % the fcn 'foreachexcitation' employs the method of images to generate a
                     % 2D field sol'n,  Inputs are the full width and length of all
                     % windings, gap length, gap mult, breadth of core window,
                     % height of core window,  number of turns, and coordinates of
                     % bottom left corner of winding cross-sections (origin at inside
                     % of core).  The  matrix of int_B2 is the output
                     
                 int_B2(i,j,:)=foreachexcitation(a,b,gaplen,gapmult, xzero,yzero,bw,h,agap,bgap,xzerogap, yzerogap, xdiv,ydiv,excwind,numwindplusgap,Ximage,Yimage);     
                  %Now the (i,j) position on all the pages has just been calculated.  Store these in the paged matrix.
                  % And define all the (j,i) positions as well.

                 int_B2(j,i,:)=int_B2(i,j,:);
         
                 excwind=zeros(numwind,1);
                 
            end
      end % end i loop
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of B2avgcalc2d%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of foreachexcitation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  function foreachexcitation changes winding excitations in sequence to obtain all field values described
 %   in the SFD method.  inputs are winding positions and locations in window, gap length and number, breadth and
 %   height of the winding window, gap locations, winding divisions for field computation, vector describing excitation and numnber
 %   of images to be computed
 % h - height of core window
% bw - breadth of winding window
% a - Height of Winding Cross Section for each winding (in h direction)
% b - Width of Winding Cross Section for each winding (in bw direction)
% xzero - Location of Center of Winding Cross Section for each winding (in h direction)
% yzero - Location of Center of Winding Cross Section for each winding (in bw direction)
% gaplen - length of gaps
% gaploc - the location of the gaps
% xdiv - the number of divisions in the x dimension for field calc
% ydiv - the number of divisions in the y dimension for field calc
% Ximage - the number of images in the x directions for method od images
% Yimage - the number of images in the Y directions for method of images
% gapmult - a gap multiplier determine in B2avgcalc2d
% xzerogap - from function ungap
% yzerogap - from function ungap
% excwind - keeps tract of which winding is excited
% numwindplusgap - keeps tract of the number of windings and the number of
% gaps for field calc
 
function [int_B2]=foreachexcitation(a,b,gaplen,gapmult, xzero,yzero,bw,h,agap,bgap,xzerogap, yzerogap, xdiv,ydiv,excwind,numwindplusgap,Ximage,Yimage)      %#ok<*INUSL>
% R. Scott Mongrain - This feels more like a debug comment.
% disp(['Now computing field for the following pattern of excitation: ' num2str(excwind')])
numwind=length(a);

numexc = sum(excwind);
% set gap elements as excited

    gaps=numwindplusgap-numwind;

for i=numwind+1:numwindplusgap %for each gap excite the ribboon with -1 current
    excwind(i)=-(numexc/gaps);
end


% combine wind vect's and gap vects

a=[a agap];
b=[b bgap];
xzero=[xzero xzerogap];
yzero=[yzero yzerogap];
   
        
         % this fcn computes the field due to a uniform rectangular current
         % distribution,  inputs are rect location, size, 
                 
   
  [int_B2]=imagecalc(a,b,xzero,yzero,bw,h,xdiv,ydiv,excwind,numwind,Ximage,Yimage);
  %%%%%%%%%%%%%%%%%%%%%%%%End for each excitation%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of imagecalc%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  Function imagecalc computes all data needed to use method of images about high perbeablilty
  %   core,  the inputs are dimension and location of each of the windings in the window, breadth and
  %   height of the window, number of divisions for windings, excited winding vector, and number of 
  %   images to be computed.  The entire system of imaged excited windings are generated.
  
   function [int_B2]=imagecalc(a,b,xzero,yzero,bw,h,xdiv,ydiv,excwind,numwind,Ximage,Yimage)

% a and b are the width and height and xzero yzero are the center points

% assign excited vectors to include winding and gap locations (eliminate unexcited windings)

% define aexc,bexc, etc

aexc=a(excwind~=0);
bexc=b(excwind~=0);
xzeroexc=xzero(excwind~=0);
yzeroexc=yzero(excwind~=0);


% We will image the geometry of excited windings to obtain a ximage by yimage system and then run Rectfldcalc, each wind/rib is on a page

apage=reshape(aexc,1,1,length(aexc));
bpage=reshape(bexc,1,1,length(bexc));
xzeropage=reshape(xzeroexc,1,1,length(xzeroexc));
yzeropage=reshape(yzeroexc,1,1,length(yzeroexc));



for n=2:(Yimage-1)/2+1                                              % Image up for y's, leave x values alone
yabove=bw*(n-1)-yzeropage(1,:,:);     % distance above winding centers to core
yzeroadd=yzeropage(1,:,:)+2.*yabove;  % %define bottom row  of matrix 
xzeroadd=xzeropage(1,:,:);  
aadd=apage(1,:,:);  
badd=bpage(1,:,:);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yzeropage=[yzeroadd; yzeropage];
xzeropage=[xzeroadd; xzeropage];
apage=[aadd; apage];
bpage=[badd; bpage];
                    
end
   
for n=2:(Yimage-1)/2+1                                              % Image down for y's, leave x values alone
    %instead of length(yzp), use correct dimesnion--added v1.4
nyzp=size(yzeropage);
n1yzp=nyzp(1);
ybelow=bw*(n-2)+yzeropage(n1yzp,:,:);     % distance above winding centers to core
yzeroadd=yzeropage(n1yzp,:,:)-2.*ybelow;  % %define bottom row  of matrix 
xzeroadd=xzeropage(n-1,:,:);  
aadd=apage(n-1,:,:);  
badd=bpage(n-1,:,:);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
yzeropage=[yzeropage; yzeroadd];
xzeropage=[xzeropage; xzeroadd];
apage=[apage; aadd];
bpage=[bpage; badd];
                    
end

% At this point we have one column of winding groups (containing the original winding group) that is yimage high,
% it will be good to define the farthest left column that we will define as column 1 and the fathest right as ximage 
% and the original as (ximage+1)/2
 
for n=2:(Ximage-1)/2+1            
    
xright=h*(n-1)-xzeropage(:,n-1,:);
xzeroadd=xzeropage(:,n-1,:)+2.*xright;  % 
yzeroadd=yzeropage(:,n-1,:);  
aadd=apage(:,n-1,:);  
badd=bpage(:,n-1,:);  

yzeropage=[yzeropage yzeroadd];             %  % Image right for x's leave y values alone
xzeropage=[xzeropage xzeroadd];                  
apage=[apage aadd];
bpage=[bpage badd];

end


for n=2:(Ximage-1)/2+1
    
xleft=h*(n-2)+xzeropage(:,1,:);
xzeroadd=xzeropage(:,1,:)-2.*xleft;  % 
yzeroadd=yzeropage(:,1,:);  
aadd=apage(:,1,:);  
badd=bpage(:,1,:);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yzeropage=[yzeroadd yzeropage];
xzeropage=[xzeroadd xzeropage];                    % Image left for x's, leave y values alone
apage=[aadd apage];
bpage=[badd bpage];
end

excwindonly = excwind(excwind~=0);

for pn=1:numwind
    % define apn,bpn, (all single values) etc for the pnth winding in the system
    apn=a(pn);
    bpn=b(pn);
    xzeropn=xzero(pn);
    yzeropn=yzero(pn);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%Compute Fields for Imaged Geometry
    % run rectfldcalc over the winding pn and obtain component fields hx and hy,  only pass through excited data and pn data
    % hx and hy are matrices over pn  that give the x and y comp of field as a fcn of pos.
  
    [int_B2(:,:,pn)] = Rectfldcalc(apage,bpage,xzeropage,yzeropage,excwindonly,apn,bpn,xzeropn,yzeropn,bw,h,xdiv,ydiv) ;  
end

%%%%%%%%%%%%%%%%%%%%%%%%End of imagecalc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of rectfldcalc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  Function Rectfldcalc takes in the imaged system of excited windings and computes the area-average squared
 %   magnetic field over the cross-section of a single transformer winding.  The inputs are the location and
 %   dimensions of the imaged system, the location and dimensions of the winding being averaged over, breadth and
 %   height of the winding window, and the number of divisions of the winding.  The single area-avergae field value
 %   is the output.
 
 function [avg_B2]=Rectfldcalc(aexc,bexc,xzeroexc,yzeroexc,excwind,apn,bpn,xzeropn,yzeropn,bw,h,xdiv,ydiv)

 %  We would like to arrive at the following for the pnth winding:
% a magnitude Hx as a function of 2 coordinates x in xspace and y in yspace in pnth wind     
% a magnitude Hy as a function of 2 coordinates x in xspace and y in yspace in pnth wind    
% a nd b are the width and height and xzero yzero are the center points

%keyboard

[i,j,k]=(size(aexc)); %#ok<RHSFN> % Find number of windings/ribbons excited in system

         excwindpages=ones(i,j,k);
          for index=1:k
                     excwindpages(:,:,index)=excwind(index)*excwindpages(:,:,index);
         end
        
halflenpn=apn./2;      %define half length and width for pn winding
halfwidpn=bpn./2;

% discretize area of pn winding 

ydivpn=ydiv;  %=round(ydiv.*bpn/bw);
xdivpn=xdiv; %round(xdiv.*apn/h); % compute number of divisions

% find where pnth winding begins and ends
xstartpn=xzeropn-halflenpn;
ystartpn=yzeropn-halfwidpn;

yfinishpn=ystartpn+bpn;
xfinishpn=xstartpn+apn;


  % Define (x,y) discretizations inside pn windings
  % a and b describe center coordintes of winding/ribbons
  % obtain discretization pn winding

  yspace=bpn/ydivpn;
  xspace=apn/xdivpn;  % compute size of divisions
     
                                                                                                  %                                   n=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                 n=2    ------------- (i,3)
% Obtain corner coordinates for all excited winds/ribs where xcorner(i,j) is the x position of corner j of winding i    l (i,2)     l
% and ycorner(i,n) is the y position of corner n of winding i                                                           l           l
                                                                                                             %          l           l   
                                                                                                             %          l           l

halfheight=aexc/2;
halfwidth=bexc/2;

halfxspace=xspace/2;
halfyspace=yspace/2;

xbotleft=xzeroexc-halfheight;
ybotleft=yzeroexc-halfwidth;        
xtopleft=xzeroexc-halfheight;
ytopleft=yzeroexc+halfwidth;   
xtopright=xzeroexc+halfheight;
ytopright=yzeroexc+halfwidth;  
xbotright=xzeroexc+halfheight;
ybotright=yzeroexc-halfwidth;     
%%%%%%%%%%%%%%%%%%%%

xstartpn=xstartpn+halfxspace;  % begin at the center of the first division
ystartpn=ystartpn+halfyspace;

 %  Iterate through field pts by going across row, then start over on row above in pn,  for every fld point get 4 r's and 4 thetas
      %  pn spans from xstartpn,ystartpn  ------->   (xfinishpn,yfinishpn)
      % The values of theta are always less than 180 deg in magnitude, example 170 deg - (-170 deg) = 20 deg, not 340
      twopi=2*pi;
  for k=1:ydivpn 
      for j=1:xdivpn    % go through a row of divisions inside pn, we compute radius angle and field one fld point at a time 
          %at the center point of the divisions 
          xsense = (xstartpn+(j-1)*xspace);
          ysense = (ystartpn+(k-1)*yspace);
          yrel = ysense - yzeroexc;
          xrel = xsense - xzeroexc;
          
          rbotleft=sqrt((xbotleft-xsense).^2+((ybotleft-ysense).^2));
          rtopleft=sqrt((xtopleft-xsense).^2+((ytopleft-ysense).^2));
          rtopright=sqrt((xtopright-xsense).^2+((ytopright-ysense).^2));
          rbotright=sqrt((xbotright-xsense).^2+((ybotright-ysense).^2));
          thetabotleft=atan2((ysense-ybotleft),(xsense-xbotleft));  
          thetatopleft=atan2((ysense-ytopleft),(xsense-xtopleft));  
          thetatopright=atan2((ysense-ytopright),(xsense-xtopright));
          thetabotright=atan2((ysense-ybotright),(xsense-xbotright));
          thetadifbrbl = mod(thetabotright-thetabotleft+pi,twopi)-pi;
          thetadiftrtl= mod(thetatopright-thetatopleft+pi,twopi)-pi;
          thetadifbltl= mod(thetabotleft-thetatopleft+pi,twopi)-pi;
          thetadifbrtr=mod(thetabotright-thetatopright+pi,twopi)-pi;
          
          % excwind contains current magnitudes of windings (which is unit excitation) and ribbons
           %numimages=size(halfwidth);
           
          %Hx=zeros(numimages(1),numimages(2),length(excwind));
          %Hy=zeros(size(Hx));

          Hx=excwindpages./(8*pi*halfwidth.*halfheight).*((yrel+halfwidth).*thetadifbrbl-(yrel-halfwidth).* thetadiftrtl+...
             (xrel+halfheight).*log(rbotleft./rtopleft)-(xrel-halfheight).*log(rbotright./rtopright));
          Hy=excwindpages./(8*pi*halfwidth.*halfheight).*((xrel+halfheight).*thetadifbltl-(xrel-halfheight).*thetadifbrtr+...
             (yrel+halfwidth).*log(rbotleft./rbotright)-(yrel-halfwidth).*log(rtopleft./rtopright));
         %end
     
          Hxtot=sum(sum(sum(Hx)));    % obtain x field at (j,k) position in pn due to entire system
          Hytot=sum(sum(sum(Hy)));
          Hxoverpn(k,j)=Hxtot;
          Hyoverpn(k,j)=Hytot;
     end
  end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All the computations are performed for all the winding groups arranged consistently with the imaged geometry,  we have  computed the
    % x and y fields at each of the divisions of pn due to all the windings in all the winding groups and stored the total resulting field values
    % in Hxoverpn and Hyoverpn,  we now have a couple of 2 dim matrix which spans pn and x and y fields at each division
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Obtain int_B2 by adding up area squared-field products 

dA=xspace*yspace;  % incremental area inside winding pn
u0 = 4e-7*pi;
Bx=u0*Hxoverpn;
By=u0*Hyoverpn;
Bsquared=(Bx.^2+By.^2) ;
%Bsquared=B.^2;  % obtain Bsquared for each increment of area inside pn 
BsquareddA=dA.*Bsquared;
int_B2=sum(sum(BsquareddA));
avg_B2=int_B2/(apn*bpn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end rectfldcalc%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%function ungap%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [agap,bgap,xzerogap,yzerogap]=ungap(gaplen,gapmult,gaploc,gapraise,h)
%this function removes gaps in a design and repaces with a ribbon of
%current of field calc.  All the inputs were calculated in B2avgcalc2d where function is called.
                
       if gapmult==1  % we have 1 center leg gap (amend a,b,xzero,yzero)
          bgap=gaplen;  
          agap=1e-6;   % width of gap is 0
          yzerogap=gapraise;
          xzerogap=0;
            
       elseif gapmult==.5;  % we have 2 gaps in the outer legs only consider the first for our window
              bgap=[gaplen gaplen];  % width of gap is 0
              yzerogap=[gapraise gapraise];
              agap=[1e-6 1e-6];
              xzerogap=[0 h];
                     
        else agap=[];
             bgap=[];
             xzerogap=[];
             yzerogap=[];
                    
                     
       end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%end ungap%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%started adding new functions here 6/23/03%%%%%
%this function plots the current waveform for the local version of litzopt

function f1 = PWLcurrentplotL(dt, numwind, I) %#ok<*STOUT>

%plots current waveform
%the input:
% I - the amplitude of the waveform at the times given in the dt vector
% dt - the vector of time segments
% numwind - the number of windings

 t = [0 cumsum(dt)];

figure('name','Current Waveforms')

for i=1:numwind   

    %do 'numwind' Current vs time plots above each other
 
    subplot(numwind,1,i), plot(t,I(i,:),'linewidth',2.5);  
    ylbl=['Winding ' num2str(i)];                          
    ylabel(ylbl, 'fontsize', 12);
    grid;
    hold on;
end


  subplot(numwind,1,numwind),xlabel('time [sec]', 'fontsize', 12);
  subplot(numwind,1,1),title('Current Waveform for Each Winding [Amps]', 'fontsize', 14);
  
  hold off;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end PWLcurrentplotL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start imscurrent%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [Ims, Irms] = Imscurrent(I, dt, numtime, numwind)

% compute mean-square current Ims and root mean squared current Irms
%the input:
% I - the amplitude of the waveform at the times given in the dt vector
% dt - the vector of time segments
% numtime - the number of time segments
% numwind - the number of windings
Ims = zeros(numwind,1);

for i = 1:numtime
    Ims = Ims + ((I(:,i+1)-I(:,i)).^2./3 + (I(:,i)+I(:,i+1)).^2)./4./sum(dt).*dt(i);
end 

% R. Scott Mongrain - Already computing RMS elsewhere, commenting display.
Irms = sqrt(Ims); %#ok<NOPRT> % Irms is obtained from the square root of the means square current

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end Imscurrent%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start find_k_ell%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function kl = find_k_ell(numwind,numtime,I,Ims,dt,N,len,rhoc,B2avg)

%calculate di/dt for each winding and approximate loss by SFD method
%the input:
% I - the amplitude of the waveform at the times given in the dt vector
% dt - the vector of time segments
% numtime - the number of time segments
% numwind - the number of windings
% Ims - the mean squared current
% N - the number of turns
% len - the average length of a turn
% rhoc - the resistivity of copper
% B2avg - the avgerage of the squared field

avgcont = zeros(numwind,numwind);
% calculate time rate of change of current for all time segments 
% and add up contibutions and store their spatial average into 'avgcont'
%dv = (I(2)-I(1))/dt(1)

for i = 1:numtime
    dv = (I(:,i+1)-I(:,i))/dt(i);
    didt2 = dv*dv';
    lm = didt2*dt(i)/sum(dt);
    avgcont = avgcont+(lm);
end

num=zeros(size(len));

% Find Modified Dynamic Loss Matrix describing each winding


for i=1:numwind
    gammahat(i)=N(i)*len(i)/pi/rhoc/4;
	Dhat(:,:,i) = gammahat(i)* B2avg(:,:,i); 
	matrix(:,:,i) = (avgcont.*Dhat(:,:,i));
	num(i)=sum(sum(matrix(:,:,i)));
end
%note: num is this number [di1/dt  di2/dt] D_jhat [di1/dt; di2/dt]

% kl is defined by  [di1/dt  di2/dt] D_jhat [di1/dt]

%                                           [di2/dt]  

%                    -------------------------------                       

%                             I_rms^2 l_w rhoc         

% we know that l_t*N=l_w
kl = num./(Ims.*len.*N*rhoc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%rnd find_k_ell%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start find optimal designs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
function [cost44, totcost,totloss,lossW,ncorrect,fp,totcost2,totloss2,loss2,n2,totcost3,totloss3,loss3,n3, awground, n2round,totcost2round, totloss2round, dia] = find_optimal_designs(awg,numwind,kl,Irms,len,N,rhoc,bb,a,insulbuild,userfp)
   %find the designs for wire sized in vector awg, and also the set of full-bobbin designs (all the "2" ones)
   %numwind, kl, Irms, len, N, rhoc are all pretty obvious.
   %userfp is the maximum packing factor relative to perfect square packing of cylinders that the user specifies; used for finding
   %full-bobbin designs.
   
%the input:
% awg - a vector of wire sizes
% kl - the value of determined in function find_k_ell
% usefp - the user specificed packing factor
% numwind - the number of windings
% Irms - the root mean squared current
% N - the number of turns
% len - the average length of a turn
% rhoc - the resistivity of copper
% insulbuild  - specifices if the insulation build on the wire is single or heavy build
% bb - breadth of bobbin window in meters
% a - a vector that conatins the height of each winding in meters


   %CRS, 11/15/02
   %JDP implement totaldesign function 1/31/03
   %JDP implemented correct version of determining optimal designs with different size wire for each winding
   %function called by this function 
    %Findtable, fprime, fprimeinverse, ff, dlookup, stranding, optpactfactor, buildabledesign, totaldesigns, costnloss2
    %3/25/03 added function dlookupinverse and totaldesignAWG
    %7/30/03 JDP - added roundtoAWGrange function to round dia matrix to the available AWG sizes given in awg vector
    
numdesigns=length(awg);

costW=zeros(numwind,numdesigns);
lossW = costW;
dia = costW;
ncorrect = costW;
awg44 = 44; % define awg44 for use i determining the kcost44, kloss44
awgfine = min(awg)-4:.1:max(awg)+4;    % Special range for figure

%determine winding with highest Irms*N product
[IrmsNproduct] = Irms.*N;
[maxindex, masterwinding] = max(IrmsNproduct);

%call to costnloss determines kcost and kloss which equal the cost and loss for #44 wire
% R. Scott Mongrain - added else case
for i=1:numwind
   if isfinite(kl(i)) ==1
      [kcost(i,:),kloss(i,:), nk(i,:)]= costnloss(Irms(i),len(i),kl(i),N(i),awg44,rhoc);
   else
       kcost(i, :) = [];
       kloss(i, :) = [];
       nk(i, :) = [];
   end
end

%determine mastercost, masterloss, and mastern for masterwinding with highest Irms*N product
%call to costnloss fills the look up table with awgfine
% R. Scott Mongrain - removed redundant comparison (isfinite is already boolean)
%                   - added else case
i=masterwinding;
   if isfinite(kl(i))
      [mastercost, masterloss, mastern]=costnloss(Irms(i),len(i),kl(i),N(i),awg,rhoc);
      [table_cost, table_loss,table_n]=costnloss(Irms(i),len(i),kl(i),N(i),awgfine,rhoc);
   else
       mastercost = [];
       masterloss = [];
       mastern = [];
       table_cost = [];
       table_loss = [];
       table_n = [];
   end

% set up emtpy matices
costW = zeros(size(mastercost));
lossW = zeros(size(masterloss));
dia = awg;

%start building the ncorrect matrix with stranding for masterwinding
ncorrect = mastern;

%find the table to look up info
table=findtable(kcost,kloss,masterwinding,table_cost,table_loss,awgfine);
% 
% %this code prints out the lookup table for debugging
% table.deriv = table.deriv';
% table.cost4deriv = table.cost4deriv';
% table.cost = table.cost';
% table.loss = table.loss';
% table.awg = table.awg';
% table.deriv(161,1) = table.deriv(160,1);
% table.cost4deriv(161,1) = table.cost4deriv(160,1);
% lookuptable = [table.deriv table.cost4deriv table.cost table.loss table.awg];

kcostmatrix=(kcost*ones(1,numdesigns));
klossmatrix=(kloss*ones(1,numdesigns));


firstmatrix = ( (kloss(masterwinding,:)./kloss).*(kcost/kcost(masterwinding,:))  )*fprime((mastercost/kcost(masterwinding,:)),table);
%test firstmatrix for NaN
[firstmatrix, kcostmatrix, klossmatrix] = testnan(firstmatrix, kcostmatrix, klossmatrix);

costW = (kcostmatrix).*(fprimeinverse(firstmatrix,table));
%test costW for NaN
[costW, kcostmatrix, klossmatrix] = testnan(costW, kcostmatrix, klossmatrix);

lossW = klossmatrix.*(ff(costW./kcostmatrix,table));
%test lossW for NaN
[lossW, costW, kcostmatrix, klossmatrix] = testnanloss(lossW, costW, kcostmatrix, klossmatrix);

dia = dlookup((costW./kcostmatrix),table);
%test costW for NaN
[dia, costW, lossW, kcostmatrix, klossmatrix] = testnandia(dia, costW, lossW, kcostmatrix, klossmatrix);

%determine ncorrect, the correct number of stranding for each windings and design
ncorrect = stranding(Irms,len,kl,N,rhoc,dia,lossW);

% round the dia matrix, which is the matrix of the wire gauge sizes for each winding to AWG sizes in range
[diaround] = roundtoAWGrange(dia, awg);

% find the cost for the rounded designs
[row,numdesignround] = size(diaround);
kcostmatrix = kcost*ones(1,numdesignround);
costWroundD = (kcostmatrix).*(dlookupinverse(diaround,table));

%find the losses for the rounded wire guage matrix
klossmatrix = kloss*ones(1, numdesignround);
lossWroundD = klossmatrix.*ff(costWroundD./kcostmatrix,table);

%calculate the stranding with the diaround, the matrix of rounded diameters and loss matrix
ncorrectALL = stranding(Irms,len,kl,N,rhoc,diaround,lossWroundD);

%determine cost44 for determing relative cost
cost44=sum(kcost);

%this determines the packing factor
fp = optpackfactor(diaround, N, ncorrectALL, bb, a, insulbuild);

%Now do the buildable designs
n2 = buildabledesign(diaround, a, ncorrectALL, fp, N, userfp, bb, insulbuild);

%the totaldesigns function rounds the number of strands in each winding
[awground, n2round] = totaldesigns(n2, diaround, numwind);

% % find the cost for the final designs with rounded awg and stranding
[cost2,loss2]= costnloss2(Irms,len,n2round,N,awground,rhoc,kl);

 %Now do the full-bobbin designs
n3 = fullbobbindesigns(awg, a, ncorrect, fp, N, userfp, bb, insulbuild);

%loop through different windings to obtain cost and loss vectors for full bobbin designs 
for i=1:numwind
    if isfinite(kl(i)) ==1
    [cost3(i,:),loss3(i,:)]= costnloss2(Irms(i),len(i),n3(i,:),N(i),awg,rhoc,kl(i));
    end
 end

 if numwind>1
   totcost=sum(costW)/cost44;                   %total cost for optimal design
   totloss=sum(lossW);                          %total loss for optimal design
   totcost2=sum(cost2)/cost44;                 %total cost for buildable design
   totloss2=sum(loss2);                         %total loss for buildable design
   totcost2round=sum(cost2)/cost44;             %total cost for ROUNDED buildable design
   totloss2round=sum(loss2);                    %total loss for ROUNDED buildable design
   totcost3=sum(cost3)/cost44;                  %total cost for buildable design
   totloss3=sum(loss3);                         %total loss for buildable design
else
    totcost=costW/cost44;           
    totloss=lossW;
    totcost2=cost2/cost44;
    totloss2=loss2;
    totcost2round=cost2/cost44;
    totloss2round=loss2;
    totcost3=cost3/cost44;
    totloss3=loss3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end find optimal designs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start find potimal design curvedata%%%%%%%%%%%%%%%%%%%%%%%%
function [totcostfine,totlossfine,totcost2fine,totloss2fine] = find_optimal_designs_curvedata(awg,numwind,kl,Irms,len,N,rhoc,bb,a,insulbuild,userfp)
   %find the designs for wire sized in vector awg, and also the set of full-bobbin designs (all the "2" ones)
   %numwind, kl, Irms, len, N, rhoc are all pretty obvious.
   %userfp is the maximum packing factor relative to perfect square packing of cylinders that the user specifies; used for finding full-bobbin designs.
  
%the inputs:
% awg - a vector of wire sizes
% kl - the value of determined in function find_k_ell
% usefp - the user specificed packing factor
% numwind - the number of windings
% Irms - the root mean squared current
% N - the number of turns
% len - the average length of a turn
% rhoc - the resistivity of copper
% insulbuild  - specifices if the insulation build on the wire is single or heavy build
% bb - breadth of bobbin window in meters
% a - a vector that conatins the height of each winding in meters

   %JDP implement totaldesign function 1/31/03
   %JDP implemented correct version of determining optimal designs with different size wire for each winding
   %function called by this function 
    %Findtable, fprime, fprimeinverse, ff, dlookup, stranding, optpactfactor, buildabledesign, totaldesigns, costnloss2
    %3/25/03 added function dlookupinverse and totaldesignAWG
    
numdesigns=length(awg);

costW=zeros(numwind,numdesigns);
lossW = costW;
dia = costW;
ncorrect = costW;
awg44 = 44; % define awg44 for use i determining the kcost44, kloss44
awgfine = min(awg):.1:max(awg);    % Special range for figure

%determine winding with highest Irms*N product
[IrmsNproduct] = Irms.*N;
[maxindex, masterwinding] = max(IrmsNproduct);

%call to costnloss determines kcost and kloss which equal the cost/loss for #44 wire
for i=1:numwind
   if isfinite(kl(i)) ==1
      [kcost(i,:),kloss(i,:), nk(i,:)]= costnloss(Irms(i),len(i),kl(i),N(i),awg44,rhoc);
   end
end

%determine mastercost, masterloss, and mastern for masterwinding with highest Irms*N product
%call to costnloss fills the look up table with awgfine
i=masterwinding;
   if isfinite(kl(i)) ==1
      [mastercost, masterloss, mastern]=costnloss(Irms(i),len(i),kl(i),N(i),awg,rhoc);
      [table_cost, table_loss, table_n]=costnloss(Irms(i),len(i),kl(i),N(i),awgfine,rhoc);
   end
   
% set up emtpy matices
costW = zeros(size(mastercost));
lossW = zeros(size(masterloss));
dia = awg;

%start building the ncorrect matrix with stranding for masterwinding
ncorrect = mastern;

%find the table to look up info
table=findtable(kcost,kloss,masterwinding,table_cost,table_loss,awgfine);

% %this code prints out the lookup table for debugging
% table.deriv = table.deriv';
% table.cost4deriv = table.cost4deriv';
% table.cost = table.cost';
% table.loss = table.loss';
% table.awg = table.awg';
% table.deriv(161,1) = table.deriv(160,1);
% table.cost4deriv(161,1) = table.cost4deriv(160,1);
% lookuptable = [table.deriv table.cost4deriv table.cost table.loss table.awg];

kcostmatrix=(kcost*ones(1,numdesigns));
klossmatrix=(kloss*ones(1,numdesigns));

%calculate the normalized cost, loss and diameter for rest of windings
firstmatrix = ( (kloss(masterwinding,:)./kloss).*(kcost/kcost(masterwinding,:))  )*fprime((mastercost/kcost(masterwinding,:)),table);
costW = (kcostmatrix).*fprimeinverse(firstmatrix,table);
lossW = klossmatrix.*ff(costW./kcostmatrix,table);
dia = dlookup(costW./kcostmatrix,table);

%determine ncorrect, the correct number of stranding for each windings and design
ncorrect = stranding(Irms,len,kl,N,rhoc,dia,lossW);

% round the dia, the matrix of the wire gauge size for each winding 
[diaround] = totaldesignsAWG(dia, awg);

% find the cost for the rounded designs
[row,numdesignround] = size(dia);
kcostmatrix = kcost*ones(1,numdesignround);
costWroundD = (kcostmatrix).*dlookupinverse(dia,table);

%find the losses for the rounded wire guage matrix
klossmatrix = kloss*ones(1, numdesignround);
lossWroundD = klossmatrix.*ff(costWroundD./kcostmatrix,table);

%calculate the stranding with the diaround, the matrix of rounded diameters and loss matrix
ncorrectALL = stranding(Irms,len,kl,N,rhoc,diaround,lossWroundD);

%determine cost44 for determing relative cost
cost44=sum(kcost);

%this determines the packing factor
fp = optpackfactor(dia, N, ncorrect, bb, a, insulbuild);

%Now do the buildable designs
n2 = buildabledesign(dia, a, ncorrect, fp, N, userfp, bb, insulbuild);

% % find the cost for the final designs with rounded awg and stranding

[cost2,loss2]= costnloss2(Irms,len,n2,N,dia,rhoc,kl);


 %Now do the full-bobbin designs
n3 = fullbobbindesigns(awg, a, ncorrect, fp, N, userfp, bb, insulbuild);

%loop through different windings to obtain cost and loss vectors for full bobbin designs 
for i=1:numwind
    if isfinite(kl(i)) ==1
    [cost3(i,:),loss3(i,:)]= costnloss2(Irms(i),len(i),n3(i,:),N(i),awg,rhoc,kl(i));
    end
 end
 
 
 if numwind>1
   totcostfine=sum(costW)/cost44;                   %total cost for optimal design
   totlossfine=sum(lossW);                          %total loss for optimal design
   totcost2fine=sum(cost2)/cost44;                  %total cost for buildable design
   totloss2fine=sum(loss2);                         %total loss for buildable design\
   totcost2round=sum(cost2)/cost44;             %total cost for ROUNDED buildable design
   totloss2round=sum(loss2);                    %total loss for ROUNDED buildable design

else
    totcostfine=costW/cost44;           
    totlossfine=lossW;
    totcost2fine=cost2/cost44;
    totloss2fine=loss2;
    totcost2roundfine=cost2/cost44;
    totloss2roundfine=loss2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end find optmal design curvedata%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start design labels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
function designlabels = designlabels(dia)

%this function creates the labels for the designs
% built by JDP on 5/28/03
%the inputs:
% dia - a vector of diameters

[row, design] = size(dia);
count = 0;
designlabels = [];

d = {'d'};

for i = 1:design
       count = count+1;
       label=num2str(count);
       templabel = strcat(d,label);
       designlabels = [designlabels; templabel];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end design labels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start plotdesignsL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f2 = plotdesignsL(tot_orig_cost, tot_orig_loss, bbdesignlabels, optdesignlabels, dia, awground, totloss, totcost, totcost2round, totloss2round, totcost3, totloss3, totcostfine,totlossfine,totcost2fine,totloss2fine) %#ok<*INUSD>
%Plot results:plots optimal design frontier
% tot_orig_cost - the relative cost for the original winding design
% tot_orig_loss - the total loss for the original winding design
% bbdesignlabels - the design labels for the buildable design
% optdesignlabels - the design labels for the optimal designs
% dia - a vector fo diameters
% awground - a vector of rounded awg sizes
% totloss - the total loss for hypotheitcal optimum design - rough data
% totcost -  the realtive cost for hypotheitcal optimum design - rough data
% totcost2round - the relative cost for the buildable designs
% totloss2round - the total loss for the  buildable designs
% totcost3 - the relative cost of the full bobbin designs
% totloss3 - the total loss for the full bobbin designs
% totcostfine - a vector of cost results for hypothetical designs with great resolution - used for the line
% totlossfine -  a vector of loss results for hypothetical designs with great resolution - used for the line
% totcost2fine -  a vector that plots the cost results for buildable designs with great resolution
% totloss2fine -  a vector that plots the loss results for buildable designs with great resolution

[row, numdesigns] = size(dia);
[rowbb, numbbdesigns] = size(awground);
figure('name','Optimal Design Frontier')

% totcostfine = totcostfine;
% totcost = totcost;
% totcost2round = totcost2round;
% tot_orig_cost = tot_orig_cost;

%use_opt = logical(totloss == totloss2);  %where the bobbin is underfilled 
% %with fine data for smooth curves
% %first plot optimal frontier in red; will be color for unbuildable designs
loglog(totcostfine,totlossfine,'-.r','linewidth',1.5)
hold on
% %Now coarse data with markers
loglog(totcost,totloss,'r^')        %plots over-full bobbin designs

%plots orignal design on optimal design frontier
if tot_orig_cost ~= 0
    loglog(tot_orig_cost, tot_orig_loss,'b*')
end

%if you want to plot the full bossin designs, use this code
% % %now full-bobbin designs in blue
% % %%%%%%%%%%%%%%%%%%%
% loglog(totcost3,totloss3,'--b','linewidth',1.5)
% % %%%%%%%%%%%%%%%%%%%%%
% loglog(totcost3,totloss3,'bs')        %plots full bobbin designs

%and finally, the buildable designs in black
%loglog(totcost2fine,totloss2fine,'-k','linewidth',1.5)
loglog(totcost2round,totloss2round,'-k','linewidth',1.5)
%Now coarse data for markers
loglog(totcost2round,totloss2round,'ko')              %plots buildable design

%label the designs on the plot
for j= 1:numbbdesigns
     text(totcost2round(j)*1.08,totloss2round(j)*1.08,bbdesignlabels(j))
     %text(totcost3(j)*1.05,totloss3(j)*1.05,num2str(dia(1,j)'),'color','b')

 end
ylabel('Loss (W)', 'fontsize', 12),xlabel('Relative Cost', 'fontsize', 12)
set(gca,'fontsize',12)
if tot_orig_cost == 0
    legend('Hypothetical optimal designs','','Best buildable designs','');
else
    legend('Hypothetical optimal designs','Original design','Best buildable designs','');
end

grid

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end plotdesignsL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start Xsectionplot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f3 = XsectionplotL(xzero, yzero, a, b, h, bw, hb, bb)

% updated by JDP 6/24/03 to fill rectangles and label better
% this function returns a plot of the winding window
% here are the definitions of h,bw,hb,bb,a,b,xzero and yzero for reference
% h - height of core window
% bw - breadth of winding window
% hb - height of bobbin window
% bb - breadth of bobbin window
% a - Height of Winding Cross Section for each winding (in h direction)
% b - Width of Winding Cross Section for each winding (in bw direction)
% xzero - Location of Center of Winding Cross Section for each winding (in h direction)
% yzero - Location of Center of Winding Cross Section for each winding (in bw direction)

figure('name','Winding Window Cross Section');

xzero = xzero.*1000;     %converts from m to mm
yzero= yzero.*1000;      %converts from m to mm
a = a.*1000;    %converts from m to mm
b = b.*1000;    %converts from m to mm

hold off
x(1,:)=[1 0 0]; %red- first winding
x(2,:)=[0 1 0]; %green- second winding
x(3,:)=[0 0 1]; %blue- third winding
x(4,:)=[1 1 0]; %yellow- fourth winding
x(5,:)=[1 0 1]; %purplish- fifth winding
xtotal = zeros(length(a),3);
%create color vector of length numwind
for g = 1:length(a)
    %find the modulus of index wrt 5
    mod5 = mod(g,5);
    %test mod5
    if mod5 == 0
        xtotal(g,:) = x(5,:);
    elseif mod5 == 4
        xtotal(g,:) = x(4,:);
    elseif mod5 == 3
        xtotal(g,:) = x(3,:);
    elseif mod5 == 2
        xtotal(g,:) = x(2,:);
    else
        xtotal(g,:) = x(1,:);
    end
end

        
%plot windings cross section and fill rectangle
  for winding=1:length(a)
      rectangle('position',[(yzero(winding)-b(winding)/2),(xzero(winding)-a(winding)/2),b(winding),a(winding)], 'EdgeColor', xtotal(winding,:))
      fillx = [(yzero(winding)-b(winding)/2),(yzero(winding)-b(winding)/2+b(winding)),(yzero(winding)-b(winding)/2+b(winding)),(yzero(winding)-b(winding)/2)];
      filly =[(xzero(winding)-a(winding)/2),(xzero(winding)-a(winding)/2),(xzero(winding)-a(winding)/2+a(winding)),(xzero(winding)-a(winding)/2+a(winding))];
      fill(fillx, filly, xtotal(winding,:))
      hold on
  end
  % R. Scott Mongrain - corrected improper window position (which has tormented
  % me for years), and added bobbin window as well.
  rectangle('Position', [0, 0, bw*1e3, h*1e3], 'LineWidth', 2)
  rectangle('Position', [(bw - bb)*1e3/2, (h - hb)*1e3, bb*1e3, hb*1e3], ...
            'Linewidth', 1.5, 'EdgeColor', ones(1, 3)/sqrt(2))
 hold off
  
 hold on

 axis equal
 xlabel('bw-- Window Breadth (mm)', 'fontsize', 14)
 ylabel('h-- Window Height (mm)', 'fontsize', 14)
 %create legend
 
 for j=1:length(a)
     title = ['Winding '   num2str(j)  '  '];
     text(yzero(j),xzero(j),title,'horizontalalignment','center')  
 end 
 
 hold off

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end XsectionplotL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start findtable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function table=findtable(kcost,kloss,masterwinding,table_cost,table_loss,awgfine)

%this function find the lookup table for determining the optimal design where each winding 
%has different size wire gauge.  It normalizes cost and loss, then find the derivative of the cost function 
% kcost - normalized cost constant
% kloss - normalized loss constant
% masterwinding - tells which winding is the mater winding
% table_cost - the cost for the table - the same size as awgfine
% table_loss - the loss for the lookup table -  the same size as awgfine
% awgfine - a vector of wire sizes in 0.1 increments

%find normqlized cost $ loss
normalcost = table_cost./kcost(masterwinding);
normalloss = table_loss./kloss(masterwinding);

%find derivative of cost function
[deriv, cost4deriv] = numderiv(normalcost, normalloss);

table.deriv = abs(deriv);
table.cost4deriv=cost4deriv;
table.cost=normalcost;
table.loss=normalloss;
table.awg = awgfine;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end findtable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555start numderiv%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [deriv, location] = numderiv(x,y)

% this function finds the numerical derivative
% x and y can be vectors
% created 2/28/03 by Jenna Pollock

deltax = diff(x);
deltay = diff(y);

deriv = abs(deltay)./deltax;

location= (x(1:(length(x)-1)) + x(2:length(x))/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end numderiv%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start fprime%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deriv=fprime(normcost,table)
    %this function looks up the derivative in the llokup table at the the value of normcost 
    %the inputs:
    %normcost - the normalized cost
    %table - the lookup table

deriv = interp1(table.cost4deriv,table.deriv,normcost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end fprime%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555start fprimeinverse%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function normalcost = fprimeinverse(deriv,table)
%this function finds the normalized cost of given a derivative and the lookup table 
    %the inputs:
    %deriv - the derviative under consideration
    %table - the lookup table
normalcost = interp1(table.deriv,table.cost4deriv,deriv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end fprimeinverse%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start dlookup%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function awg = dlookup(normcost,table)

%this function looks up the diameter cooresponding to a normalized cost in the look up table.
    %the inputs:
    %normcost - the normalized cost
    %table - the lookup table
    
awg = interp1(table.cost,table.awg,normcost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end dlookup%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start dlookupinverse%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normcost = dlookupinverse(diaround,table)

%this function looks up the normalized cost cooresponding to a diameter(awg) in the look up table.
    %the inputs:
    %diaround - the rounded diameters
    %table - the lookup table
    
normcost = interp1(table.awg,table.cost,diaround);

%%%%%%%%%%%%%%%%%%%%%%%%%%%end dlookupinverse%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start ff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normloss=ff(normcost,table)

%this function is the forward function
%given a noralized cost, it looksup the normalized loss, normloss
    %the inputs:
    %normcost - the normalized cost
    %table - the lookup table

normloss = interp1(table.cost,table.loss,normcost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end ff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start stranding%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stranding function

% Returns vectors of costs, losses, and stranding for each of the calculated
% designs where each row designates a different design.
% Inputs are 
% 'Irms' (root mean squared current), 
% 'lt' (turn length), 
% 'N' (number of turns), 
% 'awg' (a vector of strand diameters for each of the optimized designs), 
% 'rhoc' (the resistivity of copper at operating temperature)
% 'kl' (a parameter defined by time-rate of change of currents and the modified dynamic loss matrix)
% 
% The outputs: 
% 'n' is defined to be the number of strands implemented in each design.  


function [n]= stranding(Irms,lt,kl,N,rhoc,d2,lossW2)

dc = awgtomet(d2);	%convert wire guage to meters for calculation

%we need to calculate the optimal values for Fr at different guages, so do that 
%and put result into Fr[] 


Fr=1.+1./(1.-((2*COSTm(dc))./(dc.*dCOSTm(dc)))); % gives Fr 

n = numstrands(dc,kl,Fr);	%numstrands will solve equation (2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of stranding%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of totaldesignsAWG%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [diaround] = totaldesignsAWG(dia, awg)

% this function will determine the diameter for all the rounded designs
% the criterion are for numbers with a remainer < 0.2, it is rounded down.  
% If the remainder falls between 0.2 and 0.8, it is rounded both ways
% if the remainder is >0.8 the number id rounded up
    %the inputs:
    %awg - a vector of wire sizes
    %diameters - vector of wire diameters

[numwind, numdesigns] = size(dia);
diamatrix = zeros(numwind,1);
tolerance = 0.4;
max_to_round = 100;
diafinal = [];
diaout =[];
awgout=[];
numbernew = 1;
%loop thru each wire size

for i = 1:numdesigns   %wire size
  numbernew=0;
    dialow = floor(dia(:,i));
    remainder = dia(:,i) - dialow; 
    remain = remainder;
    
    %prune remain
    
        in_middle = ((remainder >= tolerance)&(remainder <= (1-tolerance)));
        low_enough_to_care = dia(:,i) < max_to_round;
        high_enough = dia(:,i) >1;
        
        to_keep = in_middle & low_enough_to_care & high_enough;
        
        remain = remain(to_keep);
        
        diasize = size(remain);
    %end of pruning
    
    %round all down, or something like that
        diamatrix = dialow + ((remainder > (1-tolerance))|(dialow<1));
        
        diafinal = [diafinal diamatrix];
         
        diasize = size(diafinal);
        numbernew = numbernew+1;
    %Done round all down step
    
    for k = 1:( length(remain))
       
      
        [maxremain, index] =  max(remain);
        diamatrix = dialow + ((remainder >= maxremain)|(dialow<1));
        diafinal = [diafinal diamatrix];
         
        diasize = size(diafinal);
        remain(index) = 0;
        numbernew = numbernew + 1;
    end 
    
    awgout = [awgout, awg(i)*ones(1,numbernew)];
    
end
diaround= diafinal;
awground = awgout;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end totaldesignAWG%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start optpackfactor%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fp] = optpackfactor(awg, N, n, bb, a, insulbuild)

%this function will determine the packing factor,fp of optimal the optimal design considering 
%the insulation build

%the inputs:
% awg - a vector of awg wires sizes
% N - a vector of the number of turns in each winding
% n - a vector of the number of strands in each design
% bb - the breadth of the bobbin window
% a - a vector containing the height of each winding
% insulbuild - the insulation build selected by user, either single or
% heavy

[numwind, numdesigns]=size(awg);
fp=zeros(numwind,numdesigns);

if insulbuild == 's'
    diameter = singlebuild(awg);
    
else
    
    diameter = heavybuild(awg);
end

for i=1:numwind
    windingarea = bb*a(i);
    
    for j=1:numdesigns
        
        fp(i,j) = (diameter(i,j)^2)*n(i,j)*N(i)/windingarea;
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5end optpackfactor%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start buildabledesigns%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n2] = buildabledesign(awg, a, n, fp, N, userfp, bb, insulbuild)

%this function will determine the a buildable design considering the
%insulation build.  It returns the number of strands, n2, that can fit in a
%buildable design for each design

%the inputs:
% awg - a vector of awg wires sizes
% n - a vector of the number of strands in each design
% fp - the packing factor required to build an optimal design
% N - a vector of the number of turns in each winding
% userfp - the packing factor specified by the user
% bb - the breadth of the bobbin window
% a - a vector containing the height of each winding
% insulbuild - the insulation build selected by user, either single or
% heavy

[numwind, numdesigns]=size(awg);

if insulbuild == 's'
    diameter = singlebuild(awg);
    
else
    
    diameter = heavybuild(awg);
end

for i=1:numwind
    
    for j=1:numdesigns
        
        if fp(i,j)>userfp
            
            n2(i,j) = userfp*bb*a(i)/((diameter(i,j)^2)*N(i));
            
        else
            n2(i,j) = n(i,j);
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end buildable designs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start totaldesigns%%%%%%%%%%%%%%%%%%%%%%%%5
function [awground, n2round] = totaldesigns(n2, awg, numwind)

% this function will determine all the rounded designs
% the criterion are for numbers with a remainer < 0.2, it is rounded down.  
% If the remainder falls between 0.2 and 0.8, it is rounded both ways
% if the remainder is >0.8 the number id rounded up

%the inputs
% n2 - a vector of the number of strands in each buildable design
% awg - a vector of wire awg sizes
% numwind - the number of windings in the component

[numwind, numdesigns] = size(awg);
n = zeros(numwind,1);
tolerance = 0.4;
max_to_round = 10;
nfinal = [];
nout =[];
awgout=[];
numbernew = 1;
%loop thru each wire size

for i = 1:numdesigns   %wire size
  numbernew=0;
    nlow = floor(n2(:,i));
    %nout = nlow;
    remainder = n2(:,i) - nlow; 
    remain = remainder;
    
    %prune remain
    
        in_middle = ((remainder >= tolerance)&(remainder <= (1-tolerance)));
        low_enough_to_care = n2(:,i) < max_to_round;
        high_enough = n2(:,i) >1;
        
        to_keep = in_middle & low_enough_to_care & high_enough;
        
        remain = remain(to_keep);
        
        nsize = size(remain);
    %end of pruning
    
    %round all down, or something like that
        n = nlow + ((remainder > (1-tolerance))|(nlow<1));
        
        nfinal = [nfinal n];
         
        nsize = size(nfinal);
        numbernew = numbernew+1;
    %Done round all down step
    
    for k = 1:( length(remain))
      
        [maxremain, index] =  max(remain);
        n = nlow + ((remainder >= maxremain)|(nlow<1));
        nfinal = [nfinal n];
         
        nsize = size(nfinal);
        remain(index) = 0;
        numbernew = numbernew + 1;
    end 
    
    awgout = [awgout, awg(:,i)*ones(1,numbernew)];    
end

n2round= nfinal;
awground = awgout;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end totaldesign%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start costnloss2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costnloss2

%this is a modified version of Costnloss to determine the cost and losses
%for buildable designs - thus "2" subscripts
%the method remains much the same, accept the n2 matrix is passed in to the function
%and it returns just the cost2 and actual_loss2 vectors.

% Returns vectors of costs, losses, and stranding for each of the calculated
% designs where each row designates a different design.
% Inputs are 
% 'Irms' (root mean squared current)
% 'lt' (turn length)
% 'N' (number of turns)
% 'awg' (a vector of strand diameters for each of the optimized designs)
% 'rhoc' (the resistivity of copper at operating temperature)
% 'kl' (a parameter defined by time-rate of change of currents and the modified dynamic loss matrix).

% The outputs 
% cost2 - are an array called 'cost' (Cost function of litz wire is assumed
% to be proportional to '(1+(k/(DC^4))+(kk))*n*lt*N' where 'k' is 1.1e-26, 
% 'kk' is 2e-9, 'n' is the number of strands in the design in question
% 'lt' is the length of a single turn, and 'N' is the number of turns.
% where, dc is strand diameter, n is strand number, lt is turn 
% length, and N is turn number)
%and
% 'actual_loss2' (the total winding power lost in each design and is equal to 'Irms^2*lt*N*rhoc*Fr/n/(pi*dc^2/4)' where all 
% parameters are as defined above and rhoc is the resistivity of copper at 
% operating temperature), and 'n' is defined to be the number of strands implemented
% in each design.  

function [cost2,actual_loss2]= costnloss2(Irms,lt,n2,N,awg,rhoc,kl)

%we need to calculate the optimal values for Fr at different guages, so do that 
%and put result into Fr[] 
%Fr=1.+1./(1.-((2*COSTm(dc))./(dc.*dCOSTm(dc)))); % gives optimal Fr but we don't want to calculate optimal
% desgins anymore.
%n = numstrands(dc,kl,Fr);	%numstrands will solve equation (2)

%now n holds a vector of the different number of strands for optimal Fr  at 
%different wire guages 
%now calculate actual_loss 
[numwind,numdesigns] = size(awg);


dc = awgtomet(awg);	%convert wire guage to meters for calculation

Fr = [];
actual_loss = [];

klm = kl*ones(1, numdesigns);

%find Fr
Fr = 1 +klm.*n2.^2.*dc.^6.*(pi/4).^3;

%calculate the loss
actual_loss2 = (((Irms.^2).*lt.*N.*rhoc)*ones(1, numdesigns)).*Fr./n2./(pi.*dc.^2/4);
if length(lt) == 2

%this will print out the dc resistnace of the first two windings.
% Rdc1 = (lt(1).*N(1).*rhoc)./(n2(1).*(pi.*dc(1).^2/4))
% Rdc2 = (lt(2).*N(2).*rhoc)./(n2(2).*(pi.*dc(2).^2/4))

end   

% now calculate cost for the cost/loss plot

Co = 0; 

cm = COSTm(dc);
    cost2 = (Co+(cm.*(dc.^2).*n2)).*( (lt.*N)*ones(1,numdesigns) );
 
%%%%%%%End of Costnloss2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of fullbobbindesigns%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n2] = fullbobbindesigns(awg, a, n, fp, N, userfp, bb, insulbuild)

%this function will determine the full-bobbin designs using the insulation build
%the inputs:
% awg - a vector of awg wires sizes
% n - a vector of the number of strands in each design
% fp - the packing factor required to build an optimal design
% N - a vector of the number of turns in each winding
% userfp - the packing factor specified by the user
% bb - the breadth of the bobbin window
% a - a vector containing the height of each winding
% insulbuild - the insulation build selected by user, either single or
% heavy

numwind = length(a);
numdesigns=length(awg);

if insulbuild == 's'
    diameter = singlebuild(awg);
    
else
    
    diameter = heavybuild(awg);
end

for i=1:numwind
    
    for j=1:numdesigns
        
        %if fp(i,j)>userfp
            
            n2(i,j) = userfp*bb*a(i)/((diameter(j)^2)*N(i));
            
        %else
        %    n2(i,j) = n(i,j);
        %end
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end fullbobbindesigns%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55start of singlebuild%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dt] = singlebuild(awg)
 
% This function determines the total diameter (in meters) of wire with single build insulation
% using the relationship determined by Charles Sullivan in "Optimal Choice for Number
% of Strands in a Litz Wire TRansformer Winding" published in IEEE Transaction on Power
% Electronics, Vol 14, No. 2, March 1999
%for wire sizes between 30-60 AWG

%the inputs:
% awg - a vector of awg wire sizes

%the output:
% dt - a vector of the total diameter of wire with insulation

alpha = 1.12;    %constant for equation of single build insulation with AWG 40 as reference
beta = .97;    %constant for equation of single build insulation with AWG 40 as reference
dr = .000079;  %the reference diameter for function of AWG 40 in meters

dc = awgtomet(awg);
dt = dr*alpha*((dc./dr).^beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end singlebuild%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start heavybuild%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dt = heavybuild(awg)
 
% This function determines the total diameter of wire with single build insulation
% using the relationship determined by Charles Sullivan in "Optimal Choice for Number
% of Strands in a Litz Wire TRansformer Winding" published in IEEE Transaction on Power
% Electronics, Vol 14, No. 2, March 1999

%the inputs:
% awg - a vector of awg wire sizes

%the output:
% dt - a vector of the total diameter of wire with insulation

alpha = 1.24;    %constant for equation of single build insulation with AWG 40 as reference
beta = .94;    %constant for equation of single build insulation with AWG 40 as reference
dr = .000079;  %the reference diameter for function of AWG 40 in meters

dc = awgtomet(awg);
dt = dr*alpha*((dc./dr).^beta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end heavybuild%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start smartround%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r] = smartround(number,sigdigits)

%this function allows a number or vector to be rounded to a specified number of significant digits

%the inputs:
% number - a vector of the numbers you wnat to round
% sigdigits - the number of significant digits you want

%the output:
% r- a vector of rounded numbers

r = [];
[row, col] = size(number);

for i = 1:row
    for j = 1:col
        leadingdigit =floor(log10(number(i,j)))+1;
    
        n = leadingdigit - sigdigits;
    
        r(i,j) = 10.^n.*round(number(i,j).*10.^-n);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end smartround%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start optdesign table%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function htmltable = optdestable(optdesignlabels, totcost, totloss, numwind, n, dia, fp)

%generate table for optimal design and their packing factors
%the inputs:
% dia - a vector of of wire diameters
% n - a vector of the number of strands in each design
% fp - the packing factor required to build an optimal design
% totcost - the total cost of the optimal designs
% totloss - the total loss of the optimal designs
% optdesignlabels - the labels for the optimal design

% R. Scott Mongrain - removed cell array on title so that quotes don't show.
fprintf('\nOptimal Design Table:\n')

totcost = smartround(totcost, 3);
totloss = smartround(totloss, 3);
totcost = totcost';
totloss = totloss';

% gaugeheading = {'Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge'; ' (W1)', ' (W2)',' (W3)',' (W4)',' (W5)',' (W6)', ' (W7)',' (W8)',' (W9)',' (W10)'};
% gaugeheading = gaugeheading(:, 1:numwind);
% heading = {'Relative','Loss','NumStrands',' NumStrands','NumStrands','NumStrands','NumStrands','NumStrands',' NumStrands','NumStrands','NumStrands','NumStrands'; 'Cost', 'in Watts', ' (W1)', ' (W2)',' (W3)',' (W4)',' (W5)',' (W6)', ' (W7)',' (W8)',' (W9)',' (W10)'};
% heading = heading(:, 1:(2+numwind));
% heading2 = {'Packing','Packing','Packing','Packing','Packing','Packing','Packing','Packing','Packing','Packing'; 'Factor (W1)','Factor (W2)','Factor (W3)','Factor (W4)','Factor (W5)','Factor (W6)','Factor (W7)','Factor (W8)','Factor (W9)','Factor (W10)'};
% heading2 = heading2(:, 1:numwind);
% 
% heading = cat(2,gaugeheading, heading);
%heading = cat(2,heading, heading2);

n = smartround(n,3);
fp = smartround(fp,3);
n = n';
dia = dia';
fp = fp';

% R. Scott Mongrain - Building table class instead of cell array.
varNames = {'Gauge', 'NumStrands', 'RelCost', 'LossInW'};
htmltable = table(dia, n, totcost, totloss, 'VariableNames', varNames);

% htmltable = [dia totcost totloss n];
% htmltable = num2cell(htmltable);
%optdesignlabels = cellstr(optdesignlabels);
%htmltable = [optdesignlabels, htmltable];

% htmltable = [heading; htmltable];
s.htmltable = htmltable; %#ok<*STRNU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end optdesigntable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start buildabletable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bdtable = buildabletable(bbdesignlabels, totcost2round, totloss2round, numwind,n2round,awground)

%generate Html table for buildable design %uses the smartround function to round cost and loss columns
%the inputs:
% awground - the rounded awg wire sizes of buildable deigns
% n2round - a vector of the number of strands in each buildable design
% totcost2round - the total cost of the buildable designs
% totloss2round - the total loss of the buildable designs
% bbdesignlabels - the labels for the buildable designs

% R. Scott Mongrain - removed cell array on title so that quotes don't show.
%create heading for the table
fprintf('\nBuildable Design Table:\n')

bdcost = totcost2round;
bdcost = smartround(bdcost,3);
bdcost = bdcost';
% gaugeheading = {'Design','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge','Gauge'; 'Number','(W1)', '(W2)','(W3)','(W4)','(W5)','(W6)', '(W7)','(W8)','(W9)','(W10)'};
% gaugeheading = gaugeheading(:, 1:1+numwind);
% bdheading = {'Relative','Loss','NumStrands',' NumStrands','NumStrands','NumStrands','NumStrands','NumStrands',' NumStrands','NumStrands','NumStrands','NumStrands'; 'Cost', 'in Watts', '(W1)', '(W2)','(W3)','(W4)','(W5)','(W6)', '(W7)','(W8)','(W9)','(W10)'};
% bdheading = bdheading(:, 1:(2+numwind));
% bdheading = cat(2,gaugeheading,bdheading);

%get data in proper format for the table
n2round = n2round';
[totloss2round] = smartround(totloss2round, 5);
totloss2round =  totloss2round';
awground = awground';
% bdtable = [awground bdcost totloss2round n2round];
% bdtable = num2cell(bdtable);
bbdesignlabels = cellstr(bbdesignlabels);
% bdtable = [bbdesignlabels, bdtable];
% bdtable = [bdheading; bdtable];

% R. Scott Mongrain - Building table class instead of cell array.
varNames = {'Design', 'Gauge', 'NumStrands', 'RelCost', 'LossInW'};
bdtable = table(bbdesignlabels, awground, n2round, bdcost, totloss2round, 'VariableNames', varNames);

s.bdtable = bdtable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end buildabletable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start roundesigns$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nrounded, awground] = rounddesigns(awg,n2, numwind) %#ok<*DEFNU>

%this program determines how to round designs and number of designs to create
%the inputs:
% numwind - the number oif windings in the component
% awg - a vector of awg sizes
% n2 - the number of strands in each  buildable design 
%the outputs:
% nrounded - the rounded vector of strandings
% awground - the vector of rounded awg sizes

numdesign = length(awg);
count = zeros(numwind,numdesign);
high = zeros(numwind,numdesign);
low = zeros(numwind,numdesign);

% determines the remainder of buildable designs then
% determines uses the remainder then rounds down if remaindern2 < 0.25
% rounds up if remaindern2>0.75 and creates both if remaindern2 is between
% 0.25 and 0.75

for k =1:numwind
    for i=1:numdesign
        if n2(k,i) <= 1
            n2(k,i) = ceil(n2(k,i));
        elseif n2(k,i) >= 10
            n2(k,i) = round(n2(k,i));
        else
            remainder = n2(k,i) - floor(n2(k,i));
            if remainder < 0.25
                n2(k,i) = floor(n2(k,i));
            elseif remainder <0.75
                count(k,i) = 1;
                low(k,i) = floor(n2(k,i));
                high(k,i) = ceil(n2(k,i));
            else
                n2(k,i) =ceil(n2(k,i));
            end
        end
    end
end

% this section of code creates the 
% rounded n matrix and awg vectors

% totaldesigns = sum(count)+numdesign;
for k = 1:numwind
    j = 1;
    for i = 1:numdesign  
        if low(k,i) == 0 && high(k,i) == 0
            awground(k,j) = awg(i);
            nrounded(k,j) = n2(k,i);
            j = j+1;
        elseif low(k,i) ~= 0 || high(k,i) ~=0
            awground(k,j) = awg(i);
            awground(k,j+1) = awg(i);
            nrounded(k,j) = low(k,i);
            nrounded(k,j+1) = high(k,i);
            j = j+2;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end%%rounddesigns%%%%%%%%%%%%%%%%%%%%

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin testnan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inputvector, kcostmatrix, klossmatrix] = testnan(inputvector, kcostmatrix, klossmatrix)

%This function prunes matrix that contains NaN that result from trying to extrapolate
%from the lookup table
%the inputs:
% inputvector -
% kcostmatirx - the cost  matrix that you want to remove Nan from
% klossmatrix - the loss matrix that you want to remove Nan from

%the outputs:
% inputvector -
% kcostmatirx - the cost  matrix without any NaN
% klossmatrix - the loss matrix without any NaN

nanvector = isnan(inputvector);

[row, columns] = size(nanvector);
 for z = 1:columns
     for z2 = 1:row
        if nanvector(z2,z) == 1    
            inputvector(:,z) = [];
            kcostmatrix(:,z) = [];
            klossmatrix(:,z) = [];           
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end testnan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin testnanloss%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inputvector, costW, kcostmatrix, klossmatrix] = testnanloss(inputvector, costW, kcostmatrix, klossmatrix)

%This function prunes matrix that contains NaN that result from trying to extrapolate
%from the lookup table
%the inputs:
% inputvector -
% kcostmatirx - the cost  matrix that you want to remove Nan from
% klossmatrix - the loss matrix that you want to remove Nan from
% costW -  the cost  matrix that you want to remove Nan from

%the outputs:
% inputvector -
% kcostmatirx - the cost matrix without any NaN
% klossmatrix - the loss matrix without any NaN
% costW -  the cost  matrix without any Nan

count = 0;
nanvector = isnan(inputvector);
[row, columns] = size(nanvector);
 for z = 1:columns
     for z2 = 1:row
         if nanvector(z2,z) == 1
             
             inputvector(:,z-count) = [];
             costW(:,z-count) = [];
             kcostmatrix(:,z-count) = [];
             klossmatrix(:,z-count) = [];
             count = count + 1;
         end
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end testnanloss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin testnandia%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inputvector, costW, lossW, kcostmatrix, klossmatrix] = testnandia(inputvector, costW, lossW, kcostmatrix, klossmatrix)

%This function prunes matrix that contains NaN that result from trying to extrapolate
%from the lookup table

%the inputs:
% inputvector -
% kcostmatirx - the cost  matrix that you want to remove Nan from
% klossmatrix - the loss matrix that you want to remove Nan from
% costW -  the cost  matrix that you want to remove Nan from
% lossW -  the loss  matrix that you want to remove Nan from

%the outputs:
% inputvector -
% kcostmatirx - the cost matrix without any NaN
% klossmatrix - the loss matrix without any NaN
% costW -  the cost  matrix matrix without any NaN
% lossW -  the loss  matrix matrix without any NaN

count = 0;
nanvector = isnan(inputvector);
[row, columns] = size(nanvector);
 for z = 1:columns
     for z2 = 1:row
         if nanvector(z2,z) == 1
             
             inputvector(:,z-count) = [];
             costW(:,z-count) = [];
             lossW(:,z-count) = [];
             kcostmatrix(:,z-count) = [];
             klossmatrix(:,z-count) = [];
             count = count + 1;
         end
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end testnandia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start roundtoAWGrange%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roundedAWG = roundtoAWGrange(dia, awg)

% this function rounded the matrix of wire diameters to the nearest value in awg
% vector passed in
% created by JDP 7/29/03
%the inputs:
% dia - a vector of diamters in meters
% awg - the vector of awg sizes you want to round to
%the output:
% roundedAWG - a vector of rounded AWG sizes that are closest to the
% diameters passed in

%determine size of dia matrix and setup roundedAWG matrix
[numwind, numdesigns] = size(dia);
roundedAWG = zeros(numwind, numdesigns);

%loop thru dia matrix and round nearest value in awg vector
for i = 1:numwind
    for j=1:numdesigns
        [minvalue, index] = min(abs(dia(i,j)-awg));
        roundedAWG(i,j) = awg(index);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of roundtoAWGrange%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start originaldesign%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tot_orig_cost, tot_orig_loss] = originaldesign(orig_numstrands, orig_wiresize, Irms, len, N, rhoc, kl, cost44)

% find cost and loss for original component design
% returns total origal cost and total original loss

%the inputs:
% Irms - a vector of the rms value of the current in each winding
% len - the average turn length for each winding
% bw - the breadth of the winding window
% N - a vetor of the number of turns in each winding
% numwind - the number of windings in the component
% rhoc - the resistivity of copper at given temperature 
% kl - a parameter defined by time-rate of change of currents and the modified dynamic loss matrix
% cost44 - the cost used to mormalize the designs to awg 44 

if orig_numstrands ~= 0
%     Irms
%     len
%     orig_numstrands
%     N
%     orig_wiresize
%     kl
    [orig_cost, orig_loss]= costnloss2(Irms, len, orig_numstrands, N, orig_wiresize, rhoc, kl);
    tot_orig_cost = sum(orig_cost)/cost44;
    tot_orig_loss = sum(orig_loss);
else
    tot_orig_cost = 0;
    tot_orig_loss = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end orignaldesign%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start windinght1d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, hb, bb, h] = windinght1d(Irms, bw, areabb, N, numwind)

%this function is designed to determine a, hb, bb and h for the 1d versions
%the inputs:
% Irms - a vector of the rms value of the current in each winding
% bw - the breadth of the winding window
% N - a vetor of the number of turns in each winding
% numwind - the number of windings in the component

%the outputs:
% a - Width of Winding Cross Section for each winding (in h direction)
% hb - the height of the bobbin
% bb -  the breadth of the bobbin
% h - the height of the core window

%this is done by : the height of the winding = ((number of turns)(Irms)/(sum of N*Irms for allwindings))*height of bobbin

%function starts here
bb = bw;
hb = areabb/bw;
h = hb;
sumNIrms= sum(N.*Irms);

%determine ht of each winding
for i=1:numwind
    
    hwinding(i) = (N(i)*Irms(i)/sumNIrms)*hb;
    
end
a = hwinding;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end windinght1d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
