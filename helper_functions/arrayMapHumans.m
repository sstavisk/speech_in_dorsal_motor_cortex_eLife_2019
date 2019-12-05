% arrayMapHumans.m
%
% Lookup table of the physical location of the channel numbers on the 
% array.
% These are from a looking down at the array perspective.
%
%
% USAGE: emap = arrayMapHumans( arrayName )
% 
% INPUT: 
%        arrayName  string specifying which array. 
%
% OUTPUT:
%        emap       10x10 matrix lookup where the location in the matrix
%                   corresponds to physical location of the channel contained
%                   in said matrix element.
%                   0 represents empty corner.
%        anatomyRotate    how many degrees to rotate counter-clockwise to orient this
%                         map to its anatomical location where up is Medial and right is
%                         Posterior (see Sergey's LFP JNE 2014 paper Figure 3).
%        anatomyTranslate how many mm to translate [x y] to make anatomically
%                          semi-accurate (at least to give an idea) placement in a figure.
%
% Created by Sergey Stavisky on 3 January 2011
% Updated 17 December 2017 with human arrays instead of monkey arrays.
% T5 map arranged by Sharlene Flesher in December 2017 based on BlackRock excel files.

function [emap, anatomyRotate, anatomyTranslate] = arrayMapHumans( arrayName )
    switch arrayName
        case 'T5_medial'   % T5's Medial (Posterior) array, implanted 2016
            emap = [ 0	 2 	 1	 3	4	6	 8	10	14	 0 ;
                     65	66	33	34	7	9	11	12	16	18 ;
                     67	68	35	36	5	17	13	23	20	22 ;
                     69	70	37	38	48	15	19	25	27	24 ;
                     71	72	39	40	42	50	54	21	29	26 ; %:. 
                     73	74	41	43	44	46	52	62	31	28 ; %:::WIRE BUNDLE::::
                     75	76	45	47	51	56	58	60	64	30 ; %:'
                     77	78	82	49	53	55	57	59	61	32 ;
                     79	80	84	86	87	89	91	94	63	95 ;
                     0	81	83	85	88	90	92	93	96	 0 
                     ];
                           
            % RELATIVE ARRAY POSITIONS
            anatomyRotate = 8.3+270; % degrees
            anatomyTranslate = [3.43 0.338]; % mm, where location between electrodes is 400um
                 
        case 'T5_lateral' %  T5's Lateral (Anterior) array, implanted 2016
            emap = [ 0	2	1	3	4	6	8	10	14	0 ;
                     65	66	33	34	7	9	11	12	16	18 ;
                     67	68	35	36	5	17	13	23	20	22 ;
                     69	70	37	38	48	15	19	25	27	24 ; %:. 
                     71	72	39	40	42	50	54	21	29	26 ; %:::WIRE BUNDLE::::
                     73	74	41	43	44	46	52	62	31	28 ; %:'
                     75	76	45	47	51	56	58	60	64	30 ;                                 
                     77	78	82	49	53	55	57	59	61	32 ;
                     79	80	84	86	87	89	91	94	63	95 ;
                     0	81	83	85	88	90	92	93	96	0 
                     ];
                     
            % RELATIVE ARRAY POSITIONS
            anatomyRotate = -10.9+270; % degrees
            anatomyTranslate = [-3.44 -0.070]; % mm, where location between electrodes is 400um

            
            % T8 anatomy coordinates are with respect to the brain scans they use in their Lancet 2017 supplement.
            % AND raw imaging data that Dr. Jaimie Henderson imported and examined in consultation with Dr.
            % Jonathan Miller.
        case 'T8_lateral'
            emap = [0  2   1  3  4  6  8 10 14  0;
                    65 66 33 34  7  9 11 12 16 18;
                    67 68 35 36  5 17 13 23 20 22;
                    69 70 37 38 48 15 19 25 27 24; %:. 
                    71 72 39 40 42 50 54 21 29 26; %:::WIRE BUNDLE::::
                    73 74 41 43 44 46 52 62 31 28; %:'
                    75 76 45 47 51 56 58 60 64 30;
                    77 78 82 49 53 55 57 59 61 32;
                    79 80 84 86 87 89 91 94 63 95; ...
                    0  81 83 85 88 90 92 93 96  0
                    ];
                
            % RELATIVE ARRAY POSITIONS
            anatomyRotate = 0-90; % degrees
            anatomyTranslate = [-4.4115 0.0442]; % mm, where location between electrodes is 400um
        
        case 'T8_medial'
            emap = [0  2   1  3  4  6  8 10 14  0;
                    65 66 33 34  7  9 11 12 16 18;
                    67 68 35 36  5 17 13 23 20 22;
                    69 70 37 38 48 15 19 25 27 24; %:. 
                    71 72 39 40 42 50 54 21 29 26; %:::WIRE BUNDLE::::
                    73 74 41 43 44 46 52 62 31 28; %:'
                    75 76 45 47 51 56 58 60 64 30;
                    77 78 82 49 53 55 57 59 61 32;
                    79 80 84 86 87 89 91 94 63 95; ...
                    0  81 83 85 88 90 92 93 96  0
                    ];     
                
           % RELATIVE ARRAY POSITIONS
           anatomyRotate = 0-90+4.89; % degrees
           anatomyTranslate = [4.3565 -0.1144]; % mm, where location between electrodes is 400um
            
        otherwise
            error( '[%s] arrayName ''%s'' is unrecognized.', mfilename, arrayName )
    end %switch arrayName
           
end % function