function s = Sum2StopLength(selCurves,norient,nAngle,ORange2)
s = zeros(nAngle,norient);
max(selCurves)
for iST= 1:size(selCurves,1)
     row= selCurves(iST,2)+1;
     col = selCurves(iST,1)+1;
     for def = -ORange2:ORange2
         col1 = col +def;
         if(col1<1)
             col1 = col1 + norient;
         end
         if( col1>norient)
             col1 = col1 - norient;
         end
         s(row,col1)= max(s(row,col1),selCurves(iST,3));
     end
end
