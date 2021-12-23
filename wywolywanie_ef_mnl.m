clear all
clc

Design.number_of_alternatives = 2;
Design.number_of_attributes = 3;

Design.FF = 0; % 1 - full factorial
Design.EF_mnl = 1; % 1 - efficient design for MNL
Design.EF_sim = 1; % number of designs drawn 
Design.EF_prec = 0.01; % precision of precision of approximation of the second order derivative (Fischer Information)
Design.EF_cs = 8; % number of choice sets
Design.EF_error = ["D"]; % D: D-error, A: A-error, B: Utility Balance, S: S-estimate


Design.Alt_name_1 = ["A"];
Design.Alt_name_2 = ["B"];


%Design.RandomCS = 12;
%Design.orthogonal_fractional_factorial = "seq"; % "seq" - sequential (orthogonality within each alternative), "sim" - simultaneous (orthogonality across alternatives)

Design.b0_1 = [0.1];
Design.b0_2 = [missing];

Design.a1_name = ["Time (min)"];
Design.a1_type = ["C"]; % C - continous, D - discrete
Design.a1_1 = ["6", "9", "12"];
Design.a1_2 = ["6", "9", "12"];
Design.b1_1 = [0.35];
Design.b1_2 = [0.4];

Design.a2_name = ["Cost (USD)"];
Design.a2_type = ["C"];
Design.a2_1 = ["5", "10", "15"];
Design.a2_2 = ["4",  "9", "14"];
Design.b2 = [0.2];

Design.a3_name = ["Comfort"];
Design.a3_type = ["D"];
Design.a3_1 = ["Low", "Average", "High"];
Design.a3_2 = ["Low", "Average", "High"];
Design.b3_1 = [NaN, 0.1, 0.2]; % Nan - zadeklarowanie poziomu bazowego
Design.b3_2 = [NaN, 0.2, 0.3];

tic
Design = DCE_DESIGN(Design);
toc