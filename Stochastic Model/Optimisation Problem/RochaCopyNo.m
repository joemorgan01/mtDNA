function MoleculeNo = RochaCopyNo (Parameters)
MoleculeNo=normrnd(Parameters(1),0.25*Parameters(1));
MoleculeNo=round(MoleculeNo,0);
end