function res = parseAmpac(jobname)
    arcfile = [jobname,'.arc'];
    arcfid = fopen(arcfile,'r');
    
    outfile = [jobname,'.out'];
    outfid = fopen(outfile,'r');
    
    res = parseAmpacFromText(arcfid, outfid);
    
    fclose(outfid);
    fclose(arcfid);
end
      
