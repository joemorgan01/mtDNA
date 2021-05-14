function MoleculeNo = MoraesCopyNo (Parameters)
MoleculeNo=0.1*Parameters(2);
MoleculeNo=normrnd(0.1*Parameters(2),2);
MoleculeNo=round(MoleculeNo,0);
end