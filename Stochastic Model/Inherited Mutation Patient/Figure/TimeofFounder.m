Y=NaN(80,1);
Y(1,1)=y(1,1);
for i=2:80
    Y(i,1)=y(i)-y(i-1);
end

Data=NaN(444,1,80);
for i=1:80
    Data(1:Y(i,1),:,i)=i.*ones(Y(i,1),1);
end
subplot(2,2,4)
histogram(Data,16,'Normalization','probability','FaceColor',[0 0 0])
ylabel('Probability','FontSize',16)
xlabel('Time of Founder Mutation (years)','FontSize',16)
axis([0 80 0 inf]);
grid on