function [FiSMao,FiSMap,FiSMam,FiSMap1,FiSMam1,FiSMap2,FiSMam2,FiSMap3,FiSMam3,FiSMap4,FiSMam4,FiSMap5,FiSMam5] = SeparatedComponents3D(...
    phaseShift,phaseShift0,FcS1aT,FcS2aT,FcS3aT,FcS4aT,FcS5aT,FcS6aT,FcS7aT,FcS8aT,FcS9aT,FcS10aT,FcS11aT,FcS12aT,FcS13aT)
phaseShift1 = phaseShift(1);
phaseShift2 = phaseShift(2);
phaseShift3 = phaseShift(3);
phaseShift4 = phaseShift(4);
phaseShift5 = phaseShift(5);
phaseShift6 = phaseShift(6);
phaseShift7 = phaseShift(7);
phaseShift8 = phaseShift(8);
phaseShift9 = phaseShift(9);
phaseShift10 = phaseShift(10);
phaseShift11 = phaseShift(11);
phaseShift12 = phaseShift(12);
MF = 1.0;
M = 0.5*[1 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0) 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0) 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0) 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0) 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0) 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0);
    1 0.5*MF*exp(-1i*phaseShift1) 0.5*MF*exp(+1i*phaseShift1) 0.5*MF*exp(-1i*2*phaseShift1) 0.5*MF*exp(+1i*2*phaseShift1) 0.5*MF*exp(-1i*3*phaseShift1) 0.5*MF*exp(+1i*3*phaseShift1) 0.5*MF*exp(-1i*4*phaseShift1) 0.5*MF*exp(+1i*4*phaseShift1) 0.5*MF*exp(-1i*5*phaseShift1) 0.5*MF*exp(+1i*5*phaseShift1) 0.5*MF*exp(-1i*6*phaseShift1) 0.5*MF*exp(+1i*6*phaseShift1);
    1 0.5*MF*exp(-1i*phaseShift2) 0.5*MF*exp(+1i*phaseShift2) 0.5*MF*exp(-1i*2*phaseShift2) 0.5*MF*exp(+1i*2*phaseShift2) 0.5*MF*exp(-1i*3*phaseShift2) 0.5*MF*exp(+1i*3*phaseShift2) 0.5*MF*exp(-1i*4*phaseShift2) 0.5*MF*exp(+1i*4*phaseShift2) 0.5*MF*exp(-1i*5*phaseShift2) 0.5*MF*exp(+1i*5*phaseShift2) 0.5*MF*exp(-1i*6*phaseShift2) 0.5*MF*exp(+1i*6*phaseShift2);
    1 0.5*MF*exp(-1i*phaseShift3) 0.5*MF*exp(+1i*phaseShift3) 0.5*MF*exp(-1i*2*phaseShift3) 0.5*MF*exp(+1i*2*phaseShift3) 0.5*MF*exp(-1i*3*phaseShift3) 0.5*MF*exp(+1i*3*phaseShift3) 0.5*MF*exp(-1i*4*phaseShift3) 0.5*MF*exp(+1i*4*phaseShift3) 0.5*MF*exp(-1i*5*phaseShift3) 0.5*MF*exp(+1i*5*phaseShift3) 0.5*MF*exp(-1i*6*phaseShift3) 0.5*MF*exp(+1i*6*phaseShift3);
    1 0.5*MF*exp(-1i*phaseShift4) 0.5*MF*exp(+1i*phaseShift4) 0.5*MF*exp(-1i*2*phaseShift4) 0.5*MF*exp(+1i*2*phaseShift4) 0.5*MF*exp(-1i*3*phaseShift4) 0.5*MF*exp(+1i*3*phaseShift4) 0.5*MF*exp(-1i*4*phaseShift4) 0.5*MF*exp(+1i*4*phaseShift4) 0.5*MF*exp(-1i*5*phaseShift4) 0.5*MF*exp(+1i*5*phaseShift4) 0.5*MF*exp(-1i*6*phaseShift4) 0.5*MF*exp(+1i*6*phaseShift4);
    1 0.5*MF*exp(-1i*phaseShift5) 0.5*MF*exp(+1i*phaseShift5) 0.5*MF*exp(-1i*2*phaseShift5) 0.5*MF*exp(+1i*2*phaseShift5) 0.5*MF*exp(-1i*3*phaseShift5) 0.5*MF*exp(+1i*3*phaseShift5) 0.5*MF*exp(-1i*4*phaseShift5) 0.5*MF*exp(+1i*4*phaseShift5) 0.5*MF*exp(-1i*5*phaseShift5) 0.5*MF*exp(+1i*5*phaseShift5) 0.5*MF*exp(-1i*6*phaseShift5) 0.5*MF*exp(+1i*6*phaseShift5);
    1 0.5*MF*exp(-1i*phaseShift6) 0.5*MF*exp(+1i*phaseShift6) 0.5*MF*exp(-1i*2*phaseShift6) 0.5*MF*exp(+1i*2*phaseShift6) 0.5*MF*exp(-1i*3*phaseShift6) 0.5*MF*exp(+1i*3*phaseShift6) 0.5*MF*exp(-1i*4*phaseShift6) 0.5*MF*exp(+1i*4*phaseShift6) 0.5*MF*exp(-1i*5*phaseShift6) 0.5*MF*exp(+1i*5*phaseShift6) 0.5*MF*exp(-1i*6*phaseShift6) 0.5*MF*exp(+1i*6*phaseShift6);
    1 0.5*MF*exp(-1i*phaseShift7) 0.5*MF*exp(+1i*phaseShift7) 0.5*MF*exp(-1i*2*phaseShift7) 0.5*MF*exp(+1i*2*phaseShift7) 0.5*MF*exp(-1i*3*phaseShift7) 0.5*MF*exp(+1i*3*phaseShift7) 0.5*MF*exp(-1i*4*phaseShift7) 0.5*MF*exp(+1i*4*phaseShift7) 0.5*MF*exp(-1i*5*phaseShift7) 0.5*MF*exp(+1i*5*phaseShift7) 0.5*MF*exp(-1i*6*phaseShift7) 0.5*MF*exp(+1i*6*phaseShift7);
    1 0.5*MF*exp(-1i*phaseShift8) 0.5*MF*exp(+1i*phaseShift8) 0.5*MF*exp(-1i*2*phaseShift8) 0.5*MF*exp(+1i*2*phaseShift8) 0.5*MF*exp(-1i*3*phaseShift8) 0.5*MF*exp(+1i*3*phaseShift8) 0.5*MF*exp(-1i*4*phaseShift8) 0.5*MF*exp(+1i*4*phaseShift8) 0.5*MF*exp(-1i*5*phaseShift8) 0.5*MF*exp(+1i*5*phaseShift8) 0.5*MF*exp(-1i*6*phaseShift8) 0.5*MF*exp(+1i*6*phaseShift8);
    1 0.5*MF*exp(-1i*phaseShift9) 0.5*MF*exp(+1i*phaseShift9) 0.5*MF*exp(-1i*2*phaseShift9) 0.5*MF*exp(+1i*2*phaseShift9) 0.5*MF*exp(-1i*3*phaseShift9) 0.5*MF*exp(+1i*3*phaseShift9) 0.5*MF*exp(-1i*4*phaseShift9) 0.5*MF*exp(+1i*4*phaseShift9) 0.5*MF*exp(-1i*5*phaseShift9) 0.5*MF*exp(+1i*5*phaseShift9) 0.5*MF*exp(-1i*6*phaseShift9) 0.5*MF*exp(+1i*6*phaseShift9);
    1 0.5*MF*exp(-1i*phaseShift10) 0.5*MF*exp(+1i*phaseShift10) 0.5*MF*exp(-1i*2*phaseShift10) 0.5*MF*exp(+1i*2*phaseShift10) 0.5*MF*exp(-1i*3*phaseShift10) 0.5*MF*exp(+1i*3*phaseShift10) 0.5*MF*exp(-1i*4*phaseShift10) 0.5*MF*exp(+1i*4*phaseShift10) 0.5*MF*exp(-1i*5*phaseShift10) 0.5*MF*exp(+1i*5*phaseShift10) 0.5*MF*exp(-1i*6*phaseShift10) 0.5*MF*exp(+1i*6*phaseShift10);
    1 0.5*MF*exp(-1i*phaseShift11) 0.5*MF*exp(+1i*phaseShift11) 0.5*MF*exp(-1i*2*phaseShift11) 0.5*MF*exp(+1i*2*phaseShift11) 0.5*MF*exp(-1i*3*phaseShift11) 0.5*MF*exp(+1i*3*phaseShift11) 0.5*MF*exp(-1i*4*phaseShift11) 0.5*MF*exp(+1i*4*phaseShift11) 0.5*MF*exp(-1i*5*phaseShift11) 0.5*MF*exp(+1i*5*phaseShift11) 0.5*MF*exp(-1i*6*phaseShift11) 0.5*MF*exp(+1i*6*phaseShift11);
    1 0.5*MF*exp(-1i*phaseShift12) 0.5*MF*exp(+1i*phaseShift12) 0.5*MF*exp(-1i*2*phaseShift12) 0.5*MF*exp(+1i*2*phaseShift12) 0.5*MF*exp(-1i*3*phaseShift12) 0.5*MF*exp(+1i*3*phaseShift12) 0.5*MF*exp(-1i*4*phaseShift12) 0.5*MF*exp(+1i*4*phaseShift12) 0.5*MF*exp(-1i*5*phaseShift12) 0.5*MF*exp(+1i*5*phaseShift12) 0.5*MF*exp(-1i*6*phaseShift12) 0.5*MF*exp(+1i*6*phaseShift12);
    ];
Minv = inv(M);

FiSMao = Minv(1,1)*FcS1aT + Minv(1,2)*FcS2aT+Minv(1,3)*FcS3aT+Minv(1,4)*FcS4aT+Minv(1,5)*FcS5aT+Minv(1,6)*FcS6aT+Minv(1,7)*FcS7aT + Minv(1,8)*FcS8aT+Minv(1,9)*FcS9aT+Minv(1,10)*FcS10aT+Minv(1,11)*FcS11aT+Minv(1,12)*FcS12aT+Minv(1,13)*FcS13aT;
FiSMap = Minv(2,1)*FcS1aT + Minv(2,2)*FcS2aT+Minv(2,3)*FcS3aT+Minv(2,4)*FcS4aT+Minv(2,5)*FcS5aT+Minv(2,6)*FcS6aT+Minv(2,7)*FcS7aT + Minv(2,8)*FcS8aT+Minv(2,9)*FcS9aT+Minv(2,10)*FcS10aT+Minv(2,11)*FcS11aT+Minv(2,12)*FcS12aT+Minv(2,13)*FcS13aT;
FiSMam = Minv(3,1)*FcS1aT + Minv(3,2)*FcS2aT+Minv(3,3)*FcS3aT+Minv(3,4)*FcS4aT+Minv(3,5)*FcS5aT+Minv(3,6)*FcS6aT+Minv(3,7)*FcS7aT + Minv(3,8)*FcS8aT+Minv(3,9)*FcS9aT+Minv(3,10)*FcS10aT+Minv(3,11)*FcS11aT+Minv(3,12)*FcS12aT+Minv(3,13)*FcS13aT;
FiSMap1 = Minv(4,1)*FcS1aT + Minv(4,2)*FcS2aT+Minv(4,3)*FcS3aT+Minv(4,4)*FcS4aT+Minv(4,5)*FcS5aT+Minv(4,6)*FcS6aT+Minv(4,7)*FcS7aT + Minv(4,8)*FcS8aT+Minv(4,9)*FcS9aT+Minv(4,10)*FcS10aT+Minv(4,11)*FcS11aT+Minv(4,12)*FcS12aT+Minv(4,13)*FcS13aT;
FiSMam1 = Minv(5,1)*FcS1aT + Minv(5,2)*FcS2aT+Minv(5,3)*FcS3aT+Minv(5,4)*FcS4aT+Minv(5,5)*FcS5aT+Minv(5,6)*FcS6aT+Minv(5,7)*FcS7aT + Minv(5,8)*FcS8aT+Minv(5,9)*FcS9aT+Minv(5,10)*FcS10aT+Minv(5,11)*FcS11aT+Minv(5,12)*FcS12aT+Minv(5,13)*FcS13aT;
FiSMap2 = Minv(6,1)*FcS1aT + Minv(6,2)*FcS2aT+Minv(6,3)*FcS3aT+Minv(6,4)*FcS4aT+Minv(6,5)*FcS5aT+Minv(6,6)*FcS6aT+Minv(6,7)*FcS7aT + Minv(6,8)*FcS8aT+Minv(6,9)*FcS9aT+Minv(6,10)*FcS10aT+Minv(6,11)*FcS11aT+Minv(6,12)*FcS12aT+Minv(6,13)*FcS13aT;
FiSMam2 = Minv(7,1)*FcS1aT + Minv(7,2)*FcS2aT+Minv(7,3)*FcS3aT+Minv(7,4)*FcS4aT+Minv(7,5)*FcS5aT+Minv(7,6)*FcS6aT+Minv(7,7)*FcS7aT + Minv(7,8)*FcS8aT+Minv(7,9)*FcS9aT+Minv(7,10)*FcS10aT+Minv(7,11)*FcS11aT+Minv(7,12)*FcS12aT+Minv(7,13)*FcS13aT;
FiSMap3 = Minv(8,1)*FcS1aT + Minv(8,2)*FcS2aT+Minv(8,3)*FcS3aT+Minv(8,4)*FcS4aT+Minv(8,5)*FcS5aT+Minv(8,6)*FcS6aT+Minv(8,7)*FcS7aT + Minv(8,8)*FcS8aT+Minv(8,9)*FcS9aT+Minv(8,10)*FcS10aT+Minv(8,11)*FcS11aT+Minv(8,12)*FcS12aT+Minv(8,13)*FcS13aT;
FiSMam3 = Minv(9,1)*FcS1aT + Minv(9,2)*FcS2aT+Minv(9,3)*FcS3aT+Minv(9,4)*FcS4aT+Minv(9,5)*FcS5aT+Minv(9,6)*FcS6aT+Minv(9,7)*FcS7aT + Minv(9,8)*FcS8aT+Minv(9,9)*FcS9aT+Minv(9,10)*FcS10aT+Minv(9,11)*FcS11aT+Minv(9,12)*FcS12aT+Minv(9,13)*FcS13aT;
FiSMap4 = Minv(10,1)*FcS1aT + Minv(10,2)*FcS2aT+Minv(10,3)*FcS3aT+Minv(10,4)*FcS4aT+Minv(10,5)*FcS5aT+Minv(10,6)*FcS6aT+Minv(10,7)*FcS7aT + Minv(10,8)*FcS8aT+Minv(10,9)*FcS9aT+Minv(10,10)*FcS10aT+Minv(10,11)*FcS11aT+Minv(10,12)*FcS12aT+Minv(10,13)*FcS13aT;
FiSMam4 = Minv(11,1)*FcS1aT + Minv(11,2)*FcS2aT+Minv(11,3)*FcS3aT+Minv(11,4)*FcS4aT+Minv(11,5)*FcS5aT+Minv(11,6)*FcS6aT+Minv(11,7)*FcS7aT + Minv(11,8)*FcS8aT+Minv(11,9)*FcS9aT+Minv(11,10)*FcS10aT+Minv(11,11)*FcS11aT+Minv(11,12)*FcS12aT+Minv(11,13)*FcS13aT;
FiSMap5 = Minv(12,1)*FcS1aT + Minv(12,2)*FcS2aT+Minv(12,3)*FcS3aT+Minv(12,4)*FcS4aT+Minv(12,5)*FcS5aT+Minv(12,6)*FcS6aT+Minv(12,7)*FcS7aT + Minv(12,8)*FcS8aT+Minv(12,9)*FcS9aT+Minv(12,10)*FcS10aT+Minv(12,11)*FcS11aT+Minv(12,12)*FcS12aT+Minv(12,13)*FcS13aT;
FiSMam5 = Minv(13,1)*FcS1aT + Minv(13,2)*FcS2aT+Minv(13,3)*FcS3aT+Minv(13,4)*FcS4aT+Minv(13,5)*FcS5aT+Minv(13,6)*FcS6aT+Minv(13,7)*FcS7aT + Minv(13,8)*FcS8aT+Minv(13,9)*FcS9aT+Minv(13,10)*FcS10aT+Minv(13,11)*FcS11aT+Minv(13,12)*FcS12aT+Minv(13,13)*FcS13aT;
end