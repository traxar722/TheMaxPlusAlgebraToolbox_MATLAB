

 %%
  Apos=[1 1 0 0 0 0 0 0;0 0 1 0 0 1 0 0;0 0 0 1 0 0 1 0;0 0 0 0 1 0 0 1;0 0 0 0 0 0 0 0];
 Aneg=[0 0 0 0 0 0 0 0;1 0 0 0 0 1 0 0;0 1 0 0 0 0 1 0;0 0 1 1 0 0 0 1;0 0 0 0 1 0 0 0];
        t1=2;
        t2=0;
        t3=1;
        t4=0;
        t5=0;
        d1=11;
        d2=12;
        d3=7;
        M0 = [0;0;0;0;0;1;1;1];
        temporizare=[t1;t2;d1+t3;d2+t4;d3+t5;d1;d2;d3]
        AA = Petrinetwork(Apos,Aneg,M0,'Mode','ptimed' ,'TimedVec',temporizare)
        sis=MaxPlusSystem(AA)

        sis2 = MaxPlusSystem()
     %%
     clear all; clc;
              Aneg = [ 
            1     0     0     0     0     0     0     0
            0     1     1     0     0     0     0     0
            0     0     0     0     0     1     0     0
            0     0     0     0     1     0     0     0
            0     0     0     0     0     0     0     1
            0     0     0     1     0     0     1     0];
        % Output Matrix (from transitions) Ao = PNTvar.output (6 x 8)
        Apos = [
            0     1     0     1     0     0     0     0
            0     0     0     0     1     0     0     0
            1     0     0     0     0     0     0     0
            0     0     0     0     0     0     1     1
            0     0     1     0     0     0     0     0
            0     0     0     0     0     1     0     0];
        % initial marking M0 = PNTvar.initial_marking
        M0 = [1 0 1 0 0 0 0 0];
        % temporizare
        tempo = [7 4 6 5 8 2 9 3];

        AA= Petrinetwork(Apos,Aneg,M0,'Mode','ptimed' ,'TimedVec',tempo)
        siss = MaxPlusSystem(AA)
          [X,Y] = sim(siss,5);

%%
clear all;clc

A_out =[1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 ] ;


A_inp = [ 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1 ];


 marking = [0 0 0 0 1 0 0 0 0 0 1 0 1 0 0 0 0];

 timing = [0 3 4 3 3 4 0 1 2 3 5 1 3 2 0 5 3];

  pnet = Petrinetwork(A_out,A_inp,marking,'Mode','ptimed' ,'TimedVec',timing);
  siss = MaxPlusSystem(pnet);

  U_init = [0 3 10 15 16 25;0 4 4 10 20 21;4 9 13 15 22 29]

  [X,Y] = sim(siss,U_init);

%%

clear all; clc;

A = maxplus([-inf 4 -inf;2 -inf -inf;-inf 3 2]);
B= maxplus([-inf -inf -inf;-inf -inf -inf;3 -inf -inf]);
C = maxplus([5 -inf 5;-inf 6 3;11 12 11]);
% verify_connected_graph(A)
% verify_connected_graph(B)
a=verify_str_connected_graph(A)
c=verify_str_connected_graph(C)
%%
a = maxplus(1);
%b = maxplus(2);

A = maxplus([2 3; 4 5])
b = maxplus([6;7])
[x,err] = solve_Axb(A, b)
 
A * maxplus([4;3])
%% solve xAxb

A = maxplus([-inf  2    5;
             -inf  -inf 3;
             -inf  -inf -inf]);

b = maxplus([4; 6; 2]);
[x,err] = solve_xAxb(A,b)


%% topology 1

% P1M1 P1M3 P2M1 P2M2 P2M3 P3M1 P3M2 M3W1 M3W2 M1W3 M1W2 M2W2 M2W3 M2W1
%  1    2    3    4     5    6   7    8    9    10   11   12   13   14
clear all; clc


A_output = [1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1];
     
     
A_input = [0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0];


all_markings = [
    0 0 0 0 0 0 0  1 0  1 0  1 0 0;  %  (8,10,12)
    0 0 0 0 0 0 0  1 0  1 0  0 1 0;  %  (8,10,13)
    0 0 0 0 0 0 0  1 0  1 0  0 0 1;  % (8,10,14)
    0 0 0 0 0 0 0  1 0  0 1  1 0 0;  %  (8,11,12)
    0 0 0 0 0 0 0  1 0  0 1  0 1 0;  %  (8,11,13)
    0 0 0 0 0 0 0  1 0  0 1  0 0 1;  %  (8,11,14) 
    0 0 0 0 0 0 0  0 1  1 0  1 0 0;  %  (9,10,12)
    0 0 0 0 0 0 0  0 1  1 0  0 1 0;  %  (9,10,13)
    0 0 0 0 0 0 0  0 1  1 0  0 0 1;  %  (9,10,14)
    0 0 0 0 0 0 0  0 1  0 1  1 0 0;  %  (9,11,12)
    0 0 0 0 0 0 0  0 1  0 1  0 1 0;  %  (9,11,13)
    0 0 0 0 0 0 0  0 1  0 1  0 0 1   %  (9,11,14)
];

timing = [3 4 1 2 3 5 3 0 0 0 0 0 0 0];

 pnet = Petrinetwork(A_output,A_input,all_markings(2,:),'Mode','ptimed' ,'TimedVec',timing);
  siss = MaxPlusSystem(pnet);
 [X,Y] = sim(siss,5);
% 
%   A0 = siss.A{1}
%    A1 = siss.A{2};
% 
%    A0star = kstar(A0);
% 
% Aa = A0star*A1
% 
% lambda = max_eig(Aa);
 R = [  3   6   0;
      1  4    13;
      2  5  0 ];

    plotMachineGantt(siss,R,X)

%%
clear all
A0 = maxplus([-inf -inf -inf -inf -inf -inf;
              4 -inf -inf -inf -inf -inf;
              -inf -inf -inf -inf -inf 2;
              -inf 8 -inf -inf -inf -inf;
              -inf -inf -inf 3 -inf -inf;
              5 -inf -inf 9 -inf -inf]);
A1 = maxplus([-inf -inf 7 -inf -inf -inf; 
               -inf -inf -inf -inf 6 -inf; 
               -inf -inf -inf -inf -inf -inf; 
               -inf -inf -inf -inf -inf -inf;
               -inf -inf -inf -inf -inf -inf 
               -inf -inf -inf -inf -inf -inf]);

Astar = kstar(A0);
A_test = Astar * A1;
lamda =max_eig(A_test)


%% topology 2
clear all
 A_output = maxplus([1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1]);
     
 A_input = maxplus([   0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0]);

all_markings = [
    0 0 0 0 0 0 0  1 0  1 0  1 0 0;  %  (8,10,12)
    0 0 0 0 0 0 0  1 0  1 0  0 1 0;  %  (8,10,13)
    0 0 0 0 0 0 0  1 0  1 0  0 0 1;  % (8,10,14)
    0 0 0 0 0 0 0  1 0  0 1  1 0 0;  %  (8,11,12)
    0 0 0 0 0 0 0  1 0  0 1  0 1 0;  %  (8,11,13)
    0 0 0 0 0 0 0  1 0  0 1  0 0 1;  %  (8,11,14) 
    0 0 0 0 0 0 0  0 1  1 0  1 0 0;  %  (9,10,12)
    0 0 0 0 0 0 0  0 1  1 0  0 1 0;  %  (9,10,13)
    0 0 0 0 0 0 0  0 1  1 0  0 0 1;  %  (9,10,14)
    0 0 0 0 0 0 0  0 1  0 1  1 0 0;  %  (9,11,12)
    0 0 0 0 0 0 0  0 1  0 1  0 1 0;  %  (9,11,13)
    0 0 0 0 0 0 0  0 1  0 1  0 0 1   %  (9,11,14)
];

timing = [3 4 1 2 3 5 3 0 0 0 0 0 0 0];

 pnet = Petrinetwork(A_output,A_input,all_markings(2,:),'Mode','ptimed' ,'TimedVec',timing);
  siss = MaxPlusSystem(pnet);
  [X,Y] = sim(siss, 5 );

 R = [  3   6   0;
      1  4    13;
      2  5  0 ];

  plotMachineGantt(siss,R,X)
