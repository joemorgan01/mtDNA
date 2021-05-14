function Y0 = InitialConditionsv2(MoleculeNo)
%% Sample a number from beta distribution, round to 2.dp (ensures intial conditions are integers), this is the mutant %
x=betarnd(18.57,10);

Y0=[MoleculeNo*(1-x) MoleculeNo*x];
Y0=round(Y0(:,:),0);

end