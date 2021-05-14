function Y0 = InitialConditionsv2(MoleculeNo)
%% Mutation % is always 0 for these patients
x=0
Y0=[MoleculeNo*(1-x) MoleculeNo*x];

Y0=round(Y0(:,:),0);

end