function CCo = Kai2Opt3D(obj,phaseShift,FcS1aT,FcS2aT,FcS3aT,FcS4aT,FcS5aT,FcS6aT,FcS7aT,FcS8aT,FcS9aT,FcS10aT,FcS11aT,FcS12aT)
% phaseShift=gpuArray(phaseShift);
phaseShift2 = phaseShift(1);
phaseShift3 = phaseShift(2);
phaseShift4 = phaseShift(3);
phaseShift5 = phaseShift(4);
phaseShift6 = phaseShift(5);
phaseShift7 = phaseShift(6);
phaseShift8 = phaseShift(7);
phaseShift9 = phaseShift(8);
phaseShift10 = phaseShift(9);
phaseShift11 = phaseShift(10);
phaseShift12 = phaseShift(11);
phaseShift13 = phaseShift(12);

phaseShift1 = 0;
phase2 = exp(-1i*phaseShift2) - exp(-1i*phaseShift1);
phase3 = exp(-1i*phaseShift3) - exp(-1i*phaseShift1);
phase4 = exp(-1i*phaseShift4) - exp(-1i*phaseShift1);
phase5 = exp(-1i*phaseShift5) - exp(-1i*phaseShift1);
phase6 = exp(-1i*phaseShift6) - exp(-1i*phaseShift1);
phase7 = exp(-1i*phaseShift7) - exp(-1i*phaseShift1);
phase8 = exp(-1i*phaseShift8) - exp(-1i*phaseShift1);
phase9 = exp(-1i*phaseShift9) - exp(-1i*phaseShift1);
phase10 = exp(-1i*phaseShift10) - exp(-1i*phaseShift1);
phase11 = exp(-1i*phaseShift11) - exp(-1i*phaseShift1);
phase12 = exp(-1i*phaseShift12) - exp(-1i*phaseShift1);
phase13 = exp(-1i*phaseShift13) - exp(-1i*phaseShift1);

phase2_2 = exp(-1i*2*phaseShift2) - exp(-1i*phaseShift1);
phase3_2 = exp(-1i*2*phaseShift3) - exp(-1i*phaseShift1);
phase4_2 = exp(-1i*2*phaseShift4) - exp(-1i*phaseShift1);
phase5_2 = exp(-1i*2*phaseShift5) - exp(-1i*phaseShift1);
phase6_2 = exp(-1i*2*phaseShift6) - exp(-1i*phaseShift1);
phase7_2 = exp(-1i*2*phaseShift7) - exp(-1i*phaseShift1);
phase8_2 = exp(-1i*2*phaseShift8) - exp(-1i*phaseShift1);
phase9_2 = exp(-1i*2*phaseShift9) - exp(-1i*phaseShift1);
phase10_2 = exp(-1i*2*phaseShift10) - exp(-1i*phaseShift1);
phase11_2 = exp(-1i*2*phaseShift11) - exp(-1i*phaseShift1);
phase12_2 = exp(-1i*2*phaseShift12) - exp(-1i*phaseShift1);
phase13_2 = exp(-1i*2*phaseShift13) - exp(-1i*phaseShift1);

phase2_3 = exp(-1i*3*phaseShift2) - exp(-1i*phaseShift1);
phase3_3 = exp(-1i*3*phaseShift3) - exp(-1i*phaseShift1);
phase4_3 = exp(-1i*3*phaseShift4) - exp(-1i*phaseShift1);
phase5_3 = exp(-1i*3*phaseShift5) - exp(-1i*phaseShift1);
phase6_3 = exp(-1i*3*phaseShift6) - exp(-1i*phaseShift1);
phase7_3 = exp(-1i*3*phaseShift7) - exp(-1i*phaseShift1);
phase8_3 = exp(-1i*3*phaseShift8) - exp(-1i*phaseShift1);
phase9_3 = exp(-1i*3*phaseShift9) - exp(-1i*phaseShift1);
phase10_3 = exp(-1i*3*phaseShift10) - exp(-1i*phaseShift1);
phase11_3 = exp(-1i*3*phaseShift11) - exp(-1i*phaseShift1);
phase12_3 = exp(-1i*3*phaseShift12) - exp(-1i*phaseShift1);
phase13_3 = exp(-1i*3*phaseShift13) - exp(-1i*phaseShift1);

phase2_4 = exp(-1i*4*phaseShift2) - exp(-1i*phaseShift1);
phase3_4 = exp(-1i*4*phaseShift3) - exp(-1i*phaseShift1);
phase4_4 = exp(-1i*4*phaseShift4) - exp(-1i*phaseShift1);
phase5_4 = exp(-1i*4*phaseShift5) - exp(-1i*phaseShift1);
phase6_4 = exp(-1i*4*phaseShift6) - exp(-1i*phaseShift1);
phase7_4 = exp(-1i*4*phaseShift7) - exp(-1i*phaseShift1);
phase8_4 = exp(-1i*4*phaseShift8) - exp(-1i*phaseShift1);
phase9_4 = exp(-1i*4*phaseShift9) - exp(-1i*phaseShift1);
phase10_4 = exp(-1i*4*phaseShift10) - exp(-1i*phaseShift1);
phase11_4 = exp(-1i*4*phaseShift11) - exp(-1i*phaseShift1);
phase12_4 = exp(-1i*4*phaseShift12) - exp(-1i*phaseShift1);
phase13_4 = exp(-1i*4*phaseShift13) - exp(-1i*phaseShift1);

phase2_5 = exp(-1i*5*phaseShift2) - exp(-1i*phaseShift1);
phase3_5 = exp(-1i*5*phaseShift3) - exp(-1i*phaseShift1);
phase4_5 = exp(-1i*5*phaseShift4) - exp(-1i*phaseShift1);
phase5_5 = exp(-1i*5*phaseShift5) - exp(-1i*phaseShift1);
phase6_5 = exp(-1i*5*phaseShift6) - exp(-1i*phaseShift1);
phase7_5 = exp(-1i*5*phaseShift7) - exp(-1i*phaseShift1);
phase8_5 = exp(-1i*5*phaseShift8) - exp(-1i*phaseShift1);
phase9_5 = exp(-1i*5*phaseShift9) - exp(-1i*phaseShift1);
phase10_5 = exp(-1i*5*phaseShift10) - exp(-1i*phaseShift1);
phase11_5 = exp(-1i*5*phaseShift11) - exp(-1i*phaseShift1);
phase12_5 = exp(-1i*5*phaseShift12) - exp(-1i*phaseShift1);
phase13_5 = exp(-1i*5*phaseShift13) - exp(-1i*phaseShift1);

phase2_6 = exp(-1i*6*phaseShift2) - exp(-1i*phaseShift1);
phase3_6 = exp(-1i*6*phaseShift3) - exp(-1i*phaseShift1);
phase4_6 = exp(-1i*6*phaseShift4) - exp(-1i*phaseShift1);
phase5_6 = exp(-1i*6*phaseShift5) - exp(-1i*phaseShift1);
phase6_6 = exp(-1i*6*phaseShift6) - exp(-1i*phaseShift1);
phase7_6 = exp(-1i*6*phaseShift7) - exp(-1i*phaseShift1);
phase8_6 = exp(-1i*6*phaseShift8) - exp(-1i*phaseShift1);
phase9_6 = exp(-1i*6*phaseShift9) - exp(-1i*phaseShift1);
phase10_6 = exp(-1i*6*phaseShift10) - exp(-1i*phaseShift1);
phase11_6 = exp(-1i*6*phaseShift11) - exp(-1i*phaseShift1);
phase12_6 = exp(-1i*6*phaseShift12) - exp(-1i*phaseShift1);
phase13_6 = exp(-1i*6*phaseShift13) - exp(-1i*phaseShift1);

% Transformation Matrix
M = [
    phase2 conj(phase2) phase2_2 conj(phase2_2) phase2_3 conj(phase2_3) phase2_4 conj(phase2_4) phase2_5 conj(phase2_5) phase2_6 conj(phase2_6)
    phase3 conj(phase3) phase3_2 conj(phase3_2) phase3_3 conj(phase3_3) phase3_4 conj(phase3_4) phase3_5 conj(phase3_5) phase3_6 conj(phase3_6)
    phase4 conj(phase4) phase4_2 conj(phase4_2) phase4_3 conj(phase4_3) phase4_4 conj(phase4_4) phase4_5 conj(phase4_5) phase4_6 conj(phase4_6)
    phase5 conj(phase5) phase5_2 conj(phase5_2) phase5_3 conj(phase5_3) phase5_4 conj(phase5_4) phase5_5 conj(phase5_5) phase5_6 conj(phase5_6)
    phase6 conj(phase6) phase6_2 conj(phase6_2) phase6_3 conj(phase6_3) phase6_4 conj(phase6_4) phase6_5 conj(phase6_5) phase6_6 conj(phase6_6)
    phase7 conj(phase7) phase7_2 conj(phase7_2) phase7_3 conj(phase7_3) phase7_4 conj(phase7_4) phase7_5 conj(phase7_5) phase7_6 conj(phase7_6)
    phase8 conj(phase8) phase8_2 conj(phase8_2) phase8_3 conj(phase8_3) phase8_4 conj(phase8_4) phase8_5 conj(phase8_5) phase8_6 conj(phase8_6)
    phase9 conj(phase9) phase9_2 conj(phase9_2) phase9_3 conj(phase9_3) phase9_4 conj(phase9_4) phase9_5 conj(phase9_5) phase9_6 conj(phase9_6)
    phase10 conj(phase10) phase10_2 conj(phase10_2) phase10_3 conj(phase10_3) phase10_4 conj(phase10_4) phase10_5 conj(phase10_5) phase10_6 conj(phase10_6)
    phase11 conj(phase11) phase11_2 conj(phase11_2) phase11_3 conj(phase11_3) phase11_4 conj(phase11_4) phase11_5 conj(phase11_5) phase11_6 conj(phase11_6)
    phase12 conj(phase12) phase12_2 conj(phase12_2) phase12_3 conj(phase12_3) phase12_4 conj(phase12_4) phase12_5 conj(phase12_5) phase12_6 conj(phase12_6)
    phase13 conj(phase13) phase13_2 conj(phase13_2) phase13_3 conj(phase13_3) phase13_4 conj(phase13_4) phase13_5 conj(phase13_5) phase13_6 conj(phase13_6)
    ];
Minv = inv(M);

FiSMap1 = Minv(1,1)*FcS1aT + Minv(1,2)*FcS2aT+Minv(1,3)*FcS3aT+Minv(1,4)*FcS4aT+Minv(1,5)*FcS5aT+Minv(1,6)*FcS6aT+Minv(1,7)*FcS7aT + Minv(1,8)*FcS8aT+Minv(1,9)*FcS9aT+Minv(1,10)*FcS10aT+Minv(1,11)*FcS11aT+Minv(1,12)*FcS12aT;
FiSMam1 = Minv(2,1)*FcS1aT + Minv(2,2)*FcS2aT+Minv(2,3)*FcS3aT+Minv(2,4)*FcS4aT+Minv(2,5)*FcS5aT+Minv(2,6)*FcS6aT+Minv(2,7)*FcS7aT + Minv(2,8)*FcS8aT+Minv(2,9)*FcS9aT+Minv(2,10)*FcS10aT+Minv(2,11)*FcS11aT+Minv(2,12)*FcS12aT;
FiSMap2 = Minv(3,1)*FcS1aT + Minv(3,2)*FcS2aT+Minv(3,3)*FcS3aT+Minv(3,4)*FcS4aT+Minv(3,5)*FcS5aT+Minv(3,6)*FcS6aT+Minv(3,7)*FcS7aT + Minv(3,8)*FcS8aT+Minv(3,9)*FcS9aT+Minv(3,10)*FcS10aT+Minv(3,11)*FcS11aT+Minv(3,12)*FcS12aT;
FiSMam2 = Minv(4,1)*FcS1aT + Minv(4,2)*FcS2aT+Minv(4,3)*FcS3aT+Minv(4,4)*FcS4aT+Minv(4,5)*FcS5aT+Minv(4,6)*FcS6aT+Minv(4,7)*FcS7aT + Minv(4,8)*FcS8aT+Minv(4,9)*FcS9aT+Minv(4,10)*FcS10aT+Minv(4,11)*FcS11aT+Minv(4,12)*FcS12aT;
FiSMap3 = Minv(5,1)*FcS1aT + Minv(5,2)*FcS2aT+Minv(5,3)*FcS3aT+Minv(5,4)*FcS4aT+Minv(5,5)*FcS5aT+Minv(5,6)*FcS6aT+Minv(5,7)*FcS7aT + Minv(5,8)*FcS8aT+Minv(5,9)*FcS9aT+Minv(5,10)*FcS10aT+Minv(5,11)*FcS11aT+Minv(5,12)*FcS12aT;
FiSMam3 = Minv(6,1)*FcS1aT + Minv(6,2)*FcS2aT+Minv(6,3)*FcS3aT+Minv(6,4)*FcS4aT+Minv(6,5)*FcS5aT+Minv(6,6)*FcS6aT+Minv(6,7)*FcS7aT + Minv(6,8)*FcS8aT+Minv(6,9)*FcS9aT+Minv(6,10)*FcS10aT+Minv(6,11)*FcS11aT+Minv(6,12)*FcS12aT;
FiSMap4 = Minv(7,1)*FcS1aT + Minv(7,2)*FcS2aT+Minv(7,3)*FcS3aT+Minv(7,4)*FcS4aT+Minv(7,5)*FcS5aT+Minv(7,6)*FcS6aT+Minv(7,7)*FcS7aT + Minv(7,8)*FcS8aT+Minv(7,9)*FcS9aT+Minv(7,10)*FcS10aT+Minv(7,11)*FcS11aT+Minv(7,12)*FcS12aT;
FiSMam4 = Minv(8,1)*FcS1aT + Minv(8,2)*FcS2aT+Minv(8,3)*FcS3aT+Minv(8,4)*FcS4aT+Minv(8,5)*FcS5aT+Minv(8,6)*FcS6aT+Minv(8,7)*FcS7aT + Minv(8,8)*FcS8aT+Minv(8,9)*FcS9aT+Minv(8,10)*FcS10aT+Minv(8,11)*FcS11aT+Minv(8,12)*FcS12aT;
FiSMap5 = Minv(9,1)*FcS1aT + Minv(9,2)*FcS2aT+Minv(9,3)*FcS3aT+Minv(9,4)*FcS4aT+Minv(9,5)*FcS5aT+Minv(9,6)*FcS6aT+Minv(9,7)*FcS7aT + Minv(9,8)*FcS8aT+Minv(9,9)*FcS9aT+Minv(9,10)*FcS10aT+Minv(9,11)*FcS11aT+Minv(9,12)*FcS12aT;
FiSMam5 = Minv(10,1)*FcS1aT + Minv(10,2)*FcS2aT+Minv(10,3)*FcS3aT+Minv(10,4)*FcS4aT+Minv(10,5)*FcS5aT+Minv(10,6)*FcS6aT+Minv(10,7)*FcS7aT + Minv(10,8)*FcS8aT+Minv(10,9)*FcS9aT+Minv(10,10)*FcS10aT+Minv(10,11)*FcS11aT+Minv(10,12)*FcS12aT;
FiSMap6 = Minv(11,1)*FcS1aT + Minv(11,2)*FcS2aT+Minv(11,3)*FcS3aT+Minv(11,4)*FcS4aT+Minv(11,5)*FcS5aT+Minv(11,6)*FcS6aT+Minv(11,7)*FcS7aT + Minv(11,8)*FcS8aT+Minv(11,9)*FcS9aT+Minv(11,10)*FcS10aT+Minv(11,11)*FcS11aT+Minv(11,12)*FcS12aT;
FiSMam6 = Minv(12,1)*FcS1aT + Minv(12,2)*FcS2aT+Minv(12,3)*FcS3aT+Minv(12,4)*FcS4aT+Minv(12,5)*FcS5aT+Minv(12,6)*FcS6aT+Minv(12,7)*FcS7aT + Minv(12,8)*FcS8aT+Minv(12,9)*FcS9aT+Minv(12,10)*FcS10aT+Minv(12,11)*FcS11aT+Minv(12,12)*FcS12aT;
CCo = abs( sum(sum( FiSMap1.*conj(FiSMam1).*obj.G )) )+...
    abs( sum(sum( FiSMap2.*conj(FiSMam2).*obj.G )) )+...
    abs( sum(sum( FiSMap3.*conj(FiSMam3).*obj.G )) )+...
     abs( sum(sum( FiSMap4.*conj(FiSMam4).*obj.G )) )+...
    abs( sum(sum( FiSMap5.*conj(FiSMam5).*obj.G )) )+...
    abs( sum(sum( FiSMap6.*conj(FiSMam6).*obj.G )) );

end
