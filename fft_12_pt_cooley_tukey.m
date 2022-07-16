%%%% 12 point FFT implemented using Cooley-Tukey algorithm %%%%

clear all;
close all;

%% Give your input here / randomize
time_dom_orig = [
                 0.7071 + 0.7071*1j;
                 0.7071 - 0.7071*1j;
                -0.7071 + 0.7071*1j; 
                -0.7071 - 0.7071*1j;
                 0.7071 + 0.7071*1j;
                 0.7071 - 0.7071*1j;
                -0.7071 + 0.7071*1j; 
                -0.7071 - 0.7071*1j;
                 0.7071 + 0.7071*1j;
                 0.7071 - 0.7071*1j;
                -0.7071 + 0.7071*1j; 
                -0.7071 - 0.7071*1j
                ];

N = 12; k1_fd = 4; n2_td = 3;

%% 4pt DFT => 0,3,6,9 ; 1,4,7,10 ; 2,5,8,11 ;
%% 3pt DFT => 0,4,8   ; 1,5,9    ; 2,6,10   ; 3,7,11  ;

%% Rearranging

time_dom_inp_1 = [time_dom_orig(1);
                  time_dom_orig(4);
                  time_dom_orig(7);
                  time_dom_orig(10)];
              
time_dom_inp_2 = [time_dom_orig(2);
                  time_dom_orig(5);
                  time_dom_orig(8);
                  time_dom_orig(11)];
              
time_dom_inp_3 =  [time_dom_orig(3);
                time_dom_orig(6);
                time_dom_orig(9);
                time_dom_orig(12)
                ];
%% 3 times 4-pt DFT               

Twiddle_0_0 = 1;
Twiddle_0_1 = 1;
Twiddle_0_2 = 1;
Twiddle_0_3 = 1;

Twiddle_1_0 = 1;
Twiddle_1_1 = -1j;
Twiddle_1_2 = -1;
Twiddle_1_3 = 1j;

Twiddle_2_0 = 1;
Twiddle_2_1 = -1;
Twiddle_2_2 = 1;
Twiddle_2_3 = -1;

Twiddle_3_0 = 1;
Twiddle_3_1 = 1j;
Twiddle_3_2 = -1;
Twiddle_3_3 = -1j;

fft_1     = [ 
              (time_dom_inp_1(1)*Twiddle_0_0+time_dom_inp_1(2)*Twiddle_0_1+time_dom_inp_1(3)*Twiddle_0_2+time_dom_inp_1(4)*Twiddle_0_3);
              (time_dom_inp_1(1)*Twiddle_1_0+time_dom_inp_1(2)*Twiddle_1_1+time_dom_inp_1(3)*Twiddle_1_2+time_dom_inp_1(4)*Twiddle_1_3);
              (time_dom_inp_1(1)*Twiddle_2_0+time_dom_inp_1(2)*Twiddle_2_1+time_dom_inp_1(3)*Twiddle_2_2+time_dom_inp_1(4)*Twiddle_2_3);
              (time_dom_inp_1(1)*Twiddle_3_0+time_dom_inp_1(2)*Twiddle_3_1+time_dom_inp_1(3)*Twiddle_3_2+time_dom_inp_1(4)*Twiddle_3_3)
             ];

fft_2     = [ 
              (time_dom_inp_2(1)*Twiddle_0_0+time_dom_inp_2(2)*Twiddle_0_1+time_dom_inp_2(3)*Twiddle_0_2+time_dom_inp_2(4)*Twiddle_0_3);
              (time_dom_inp_2(1)*Twiddle_1_0+time_dom_inp_2(2)*Twiddle_1_1+time_dom_inp_2(3)*Twiddle_1_2+time_dom_inp_2(4)*Twiddle_1_3);
              (time_dom_inp_2(1)*Twiddle_2_0+time_dom_inp_2(2)*Twiddle_2_1+time_dom_inp_2(3)*Twiddle_2_2+time_dom_inp_2(4)*Twiddle_2_3);
              (time_dom_inp_2(1)*Twiddle_3_0+time_dom_inp_2(2)*Twiddle_3_1+time_dom_inp_2(3)*Twiddle_3_2+time_dom_inp_2(4)*Twiddle_3_3)
             ];

fft_3     = [ 
              (time_dom_inp_3(1)*Twiddle_0_0+time_dom_inp_3(2)*Twiddle_0_1+time_dom_inp_3(3)*Twiddle_0_2+time_dom_inp_3(4)*Twiddle_0_3);
              (time_dom_inp_3(1)*Twiddle_1_0+time_dom_inp_3(2)*Twiddle_1_1+time_dom_inp_3(3)*Twiddle_1_2+time_dom_inp_3(4)*Twiddle_1_3);
              (time_dom_inp_3(1)*Twiddle_2_0+time_dom_inp_3(2)*Twiddle_2_1+time_dom_inp_3(3)*Twiddle_2_2+time_dom_inp_3(4)*Twiddle_2_3);
              (time_dom_inp_3(1)*Twiddle_3_0+time_dom_inp_3(2)*Twiddle_3_1+time_dom_inp_3(3)*Twiddle_3_2+time_dom_inp_3(4)*Twiddle_3_3)
             ];

%% Twiddle Factor Multiplication       %% WN^n2,k;    n2=0..n2_td-1,k1=0..k1_td-1
twiddle_1  = 1;
twiddle_2  = 1;
twiddle_3  = 1;
twiddle_4  = 1;
twiddle_5  = exp(-pi*1j/6);
twiddle_6  = exp(-pi*1j/3);
twiddle_7  = 1;
twiddle_8  = exp(-pi*1j/3);
twiddle_9  = exp(-2*pi*1j/3);
twiddle_10 = 1;
twiddle_11 = exp(-pi*1j/2);
twiddle_12 = exp(-pi*1j);

product_array = [fft_1(1)*twiddle_1;
                 fft_2(1)*twiddle_2;
                 fft_3(1)*twiddle_3;
                 fft_1(2)*twiddle_4;
                 fft_2(2)*twiddle_5;
                 fft_3(2)*twiddle_6;
                 fft_1(3)*twiddle_7;
                 fft_2(3)*twiddle_8;
                 fft_3(3)*twiddle_9;
                 fft_1(4)*twiddle_10;
                 fft_2(4)*twiddle_11;
                 fft_3(4)*twiddle_12
                 ];
%% 4 times 3-pt DFT             

Twiddle3_0_0 = 1;
Twiddle3_0_1 = 1;
Twiddle3_0_2 = 1;

Twiddle3_1_0 = 1;
Twiddle3_1_1 = -0.5-0.866*1j;
Twiddle3_1_2 = -0.5+0.866*1j;

Twiddle3_2_0 = 1;
Twiddle3_2_1 = -0.5+0.866*1j;
Twiddle3_2_2 = -0.5-0.866*1j;

fft_2_out_1 = [ 
              (product_array(1)*Twiddle3_0_0+product_array(2)*Twiddle3_0_1+product_array(3)*Twiddle3_0_2);
              (product_array(1)*Twiddle3_1_0+product_array(2)*Twiddle3_1_1+product_array(3)*Twiddle3_1_2);
              (product_array(1)*Twiddle3_2_0+product_array(2)*Twiddle3_2_1+product_array(3)*Twiddle3_2_2)
             ];
fft_2_out_2 = [ 
              (product_array(4)*Twiddle3_0_0+product_array(5)*Twiddle3_0_1+product_array(6)*Twiddle3_0_2);
              (product_array(4)*Twiddle3_1_0+product_array(5)*Twiddle3_1_1+product_array(6)*Twiddle3_1_2);
              (product_array(4)*Twiddle3_2_0+product_array(5)*Twiddle3_2_1+product_array(6)*Twiddle3_2_2)
             ];
fft_2_out_3 = [ 
              (product_array(7)*Twiddle3_0_0+product_array(8)*Twiddle3_0_1+product_array(9)*Twiddle3_0_2);
              (product_array(7)*Twiddle3_1_0+product_array(8)*Twiddle3_1_1+product_array(9)*Twiddle3_1_2);
              (product_array(7)*Twiddle3_2_0+product_array(8)*Twiddle3_2_1+product_array(9)*Twiddle3_2_2)
             ];
fft_2_out_4 = [ 
              (product_array(10)*Twiddle3_0_0+product_array(11)*Twiddle3_0_1+product_array(12)*Twiddle3_0_2);
              (product_array(10)*Twiddle3_1_0+product_array(11)*Twiddle3_1_1+product_array(12)*Twiddle3_1_2);
              (product_array(10)*Twiddle3_2_0+product_array(11)*Twiddle3_2_1+product_array(12)*Twiddle3_2_2)
             ];

%% Output

fft_final_out = [
                 fft_2_out_1(1);
                 fft_2_out_2(1);
                 fft_2_out_3(1);
                 fft_2_out_4(1);
                 fft_2_out_1(2);
                 fft_2_out_2(2);
                 fft_2_out_3(2);
                 fft_2_out_4(2);
                 fft_2_out_1(3);
                 fft_2_out_2(3);
                 fft_2_out_3(3);
                 fft_2_out_4(3)
                 ];

%% IFFT

%% 3pt DFT => 0,4,8   ; 1,5,9    ; 2,6,10   ; 3,7,11  ;
%% 4pt DFT => 0,3,6,9 ; 1,4,7,10 ; 2,5,8,11 ;

%% Rearranging

freq_dom_rearr_0 = [
                     fft_final_out(1);
                     fft_final_out(5);
                     fft_final_out(9)
                   ];

freq_dom_rearr_1 = [
                     fft_final_out(2);
                     fft_final_out(6);
                     fft_final_out(10)
                   ];
               
freq_dom_rearr_2 = [
                     fft_final_out(3);
                     fft_final_out(7);
                     fft_final_out(11)
                   ];
               
freq_dom_rearr_3 = [
                     fft_final_out(4);
                     fft_final_out(8);
                     fft_final_out(12)
                   ];
               
%% 4-times 3pt IDFT

Twiddle3_0_0 = 1;
Twiddle3_0_1 = 1;
Twiddle3_0_2 = 1;

Twiddle3_1_0 = 1;
Twiddle3_1_1 = -0.5+0.866*1j;
Twiddle3_1_2 = -0.5-0.866*1j;

Twiddle3_2_0 = 1;
Twiddle3_2_1 = -0.5-0.866*1j;
Twiddle3_2_2 = -0.5+0.866*1j;

fft_2_out_1 = [ 
              (freq_dom_rearr_0(1)*Twiddle3_0_0+freq_dom_rearr_0(2)*Twiddle3_0_1+freq_dom_rearr_0(3)*Twiddle3_0_2);
              (freq_dom_rearr_0(1)*Twiddle3_1_0+freq_dom_rearr_0(2)*Twiddle3_1_1+freq_dom_rearr_0(3)*Twiddle3_1_2);
              (freq_dom_rearr_0(1)*Twiddle3_2_0+freq_dom_rearr_0(2)*Twiddle3_2_1+freq_dom_rearr_0(3)*Twiddle3_2_2)
             ];
fft_2_out_2 = [ 
              (freq_dom_rearr_1(1)*Twiddle3_0_0+freq_dom_rearr_1(2)*Twiddle3_0_1+freq_dom_rearr_1(3)*Twiddle3_0_2);
              (freq_dom_rearr_1(1)*Twiddle3_1_0+freq_dom_rearr_1(2)*Twiddle3_1_1+freq_dom_rearr_1(3)*Twiddle3_1_2);
              (freq_dom_rearr_1(1)*Twiddle3_2_0+freq_dom_rearr_1(2)*Twiddle3_2_1+freq_dom_rearr_1(3)*Twiddle3_2_2)
             ];
fft_2_out_3 = [ 
              (freq_dom_rearr_2(1)*Twiddle3_0_0+freq_dom_rearr_2(2)*Twiddle3_0_1+freq_dom_rearr_2(3)*Twiddle3_0_2);
              (freq_dom_rearr_2(1)*Twiddle3_1_0+freq_dom_rearr_2(2)*Twiddle3_1_1+freq_dom_rearr_2(3)*Twiddle3_1_2);
              (freq_dom_rearr_2(1)*Twiddle3_2_0+freq_dom_rearr_2(2)*Twiddle3_2_1+freq_dom_rearr_2(3)*Twiddle3_2_2)
             ];
fft_2_out_4 = [ 
              (freq_dom_rearr_3(1)*Twiddle3_0_0+freq_dom_rearr_3(2)*Twiddle3_0_1+freq_dom_rearr_3(3)*Twiddle3_0_2);
              (freq_dom_rearr_3(1)*Twiddle3_1_0+freq_dom_rearr_3(2)*Twiddle3_1_1+freq_dom_rearr_3(3)*Twiddle3_1_2);
              (freq_dom_rearr_3(1)*Twiddle3_2_0+freq_dom_rearr_3(2)*Twiddle3_2_1+freq_dom_rearr_3(3)*Twiddle3_2_2)
             ];
         
%% Twiddle Factor Multiplication       
twiddle2_1  = 1;
twiddle2_2  = 1;
twiddle2_3  = 1;
twiddle2_4  = 1;

twiddle2_5  = 1;
twiddle2_6  = 0.8660 + 0.5000j;
twiddle2_7  = 0.5000 + 0.8660j;
twiddle2_8  = 1j;

twiddle2_9  = 1;
twiddle2_10 = 0.5000 + 0.8660j;
twiddle2_11 = -0.5000 + 0.8660j;
twiddle2_12 = -1;

product2_array = [
                 fft_2_out_1(1)*twiddle2_1;
                 fft_2_out_2(1)*twiddle2_2;
                 fft_2_out_3(1)*twiddle2_3;
                 fft_2_out_4(1)*twiddle2_4;
                 fft_2_out_1(2)*twiddle2_5;
                 fft_2_out_2(2)*twiddle2_6;
                 fft_2_out_3(2)*twiddle2_7;
                 fft_2_out_4(2)*twiddle2_8;
                 fft_2_out_1(3)*twiddle2_9;
                 fft_2_out_2(3)*twiddle2_10;
                 fft_2_out_3(3)*twiddle2_11;
                 fft_2_out_4(3)*twiddle2_12
                 ];

%%
%% 3 times 4-pt DFT               

Twiddle2_0_0 = 1;
Twiddle2_0_1 = 1;
Twiddle2_0_2 = 1;
Twiddle2_0_3 = 1;

Twiddle2_1_0 = 1;
Twiddle2_1_1 = 1j;
Twiddle2_1_2 = -1;
Twiddle2_1_3 = -1j;

Twiddle2_2_0 = 1;
Twiddle2_2_1 = -1;
Twiddle2_2_2 = 1;
Twiddle2_2_3 = -1;

Twiddle2_3_0 = 1;
Twiddle2_3_1 = -1j;
Twiddle2_3_2 = -1;
Twiddle2_3_3 = 1j;

ifft_1     = [ 
              (product2_array(1)*Twiddle2_0_0+product2_array(2)*Twiddle2_0_1+product2_array(3)*Twiddle2_0_2+product2_array(4)*Twiddle2_0_3);
              (product2_array(1)*Twiddle2_1_0+product2_array(2)*Twiddle2_1_1+product2_array(3)*Twiddle2_1_2+product2_array(4)*Twiddle2_1_3);
              (product2_array(1)*Twiddle2_2_0+product2_array(2)*Twiddle2_2_1+product2_array(3)*Twiddle2_2_2+product2_array(4)*Twiddle2_2_3);
              (product2_array(1)*Twiddle2_3_0+product2_array(2)*Twiddle2_3_1+product2_array(3)*Twiddle2_3_2+product2_array(4)*Twiddle2_3_3)
             ];

ifft_2     = [ 
              (product2_array(5)*Twiddle2_0_0+product2_array(6)*Twiddle2_0_1+product2_array(7)*Twiddle2_0_2+product2_array(8)*Twiddle2_0_3);
              (product2_array(5)*Twiddle2_1_0+product2_array(6)*Twiddle2_1_1+product2_array(7)*Twiddle2_1_2+product2_array(8)*Twiddle2_1_3);
              (product2_array(5)*Twiddle2_2_0+product2_array(6)*Twiddle2_2_1+product2_array(7)*Twiddle2_2_2+product2_array(8)*Twiddle2_2_3);
              (product2_array(5)*Twiddle2_3_0+product2_array(6)*Twiddle2_3_1+product2_array(7)*Twiddle2_3_2+product2_array(8)*Twiddle2_3_3)
             ];

ifft_3     = [ 
              (product2_array(9)*Twiddle2_0_0+product2_array(10)*Twiddle2_0_1+product2_array(11)*Twiddle2_0_2+product2_array(12)*Twiddle2_0_3);
              (product2_array(9)*Twiddle2_1_0+product2_array(10)*Twiddle2_1_1+product2_array(11)*Twiddle2_1_2+product2_array(12)*Twiddle2_1_3);
              (product2_array(9)*Twiddle2_2_0+product2_array(10)*Twiddle2_2_1+product2_array(11)*Twiddle2_2_2+product2_array(12)*Twiddle2_2_3);
              (product2_array(9)*Twiddle2_3_0+product2_array(10)*Twiddle2_3_1+product2_array(11)*Twiddle2_3_2+product2_array(12)*Twiddle2_3_3)
             ];
         
%% IFFT Output

ifft_final_out = [
                 ifft_1(1);
                 ifft_2(1);
                 ifft_3(1);
                 ifft_1(2);
                 ifft_2(2);
                 ifft_3(2);
                 ifft_1(3);
                 ifft_2(3);
                 ifft_3(3);
                 ifft_1(4);
                 ifft_2(4);
                 ifft_3(4)
                 ];
             
%% Scaling

ifft_final_out = ifft_final_out/12;
