function [ medDet,polyDet ] = medDetPdcData( flux,filterLen,quarters, order)
%Median detrend pdc data that has been stitched together.
%Remove order order polynomial to determine the 
%Then do a median detrending across the entire light curve.

nquarters=unique(quarters);
index=1:length(flux);

polyDet=zeros(length(flux),1);

%figure
for nq=1:length(nquarters)
   clear short
   thisQ=nquarters(nq);
   short=flux(quarters==thisQ);
   shortindex=index(quarters==thisQ);
   c=1:length(short);
   
   pfit=polyfit(c',short,order);
   pval=polyval(pfit,c');
   
   polyDet(shortindex)=short./pval-1;
    
  % plot(short,'ro')
  % hold on
  % plot(pval,'b-')
  % hold off
  % drawnow
end


%Do the median detrending now.

medFilt=medfilt1(polyDet,filterLen); 
medDet=polyDet-medFilt;

%Test if median detrending gives out lots of the same value
%Then try a different filter
length(medDet(medDet==0))
if length(medDet(medDet==0)) > 0.02*length(medDet)
   c1=timeseries(polyDet,1:length(polyDet));
   interval=[0, (1/3)/filterLen]
   meanpolyDet=mean(polyDet);
    
   idealfilt = idealfilter(c1,interval,'pass');
   
   notchFilt=idealfilt.Data+meanpolyDet;
   medDet=polyDet-notchFilt;
end


%Plot for testing

%figure
%subplot(2,1,1)
%plot(flux,'k.-')
%subplot(2,1,2)
%plot(polyDet,'ro')
%xlim([1e4,2e4])
%hold on
%plot(medDet,'bx')
%hold off
%drawnow
end

