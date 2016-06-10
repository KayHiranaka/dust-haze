import sys
sys.path.append('/Users/ooookay/Documents/ASTRO/codes/kay-repo/Python/module/')
#sys.path.append('/Users/kayhiranaka/Documents/Kaystuff/python/modules/')
import BDdb
import asciitable
import cPickle

db = BDdb.get_db('/Users/ooookay/Dropbox/BDNYCdb/BDNYC.db')

#L6 field red ones
result = db.query.execute("SELECT a.source_id, d.designation, b.id, b.bibtex, b.shortname, c.gravity, c.spectral_type, (c1.magnitude-c2.magnitude) AS color FROM spectra a JOIN photometry c1 ON a.source_id=c1.source_id JOIN photometry c2 ON a.source_id=c2.source_id JOIN publications b ON a.publication_id=b.id JOIN spectral_types c ON a.source_id=c.source_id JOIN sources d ON a.source_id=d.id WHERE c.gravity IS NULL AND c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND a.instrument_id=6 AND c.regime='OPT' AND c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND a.instrument_id=6 AND (c.spectral_type<=15.5 AND c.spectral_type>=10.0) AND ((color>=1.33 AND (c.spectral_type=10 OR c.spectral_type=10.5)) OR (color>=1.33 AND (c.spectral_type=11 OR c.spectral_type=11.5)) OR (color>=1.67 AND c.spectral_type=12) OR (color>=1.63 AND (c.spectral_type=13 OR c.spectral_type=13.5)) OR (color>=1.86 AND (c.spectral_type=14 OR c.spectral_type=14.5)) OR (color>=1.46 AND (c.spectral_type=15 OR c.spectral_type=15.5)) OR (color>=1.89 AND (c.spectral_type=16 OR c.spectral_type=16.5)))").fetchall()
out = open('Files/dbquery_fg_L6.txt', 'wb')
cPickle.dump(result, out)
quit()

#spex prism ref for red
result = db.query.execute("SELECT a.source_id, d.designation, b.id, b.bibtex, b.shortname, c.gravity, c.spectral_type, (c1.magnitude-c2.magnitude) AS color FROM spectra a JOIN photometry c1 ON a.source_id=c1.source_id JOIN photometry c2 ON a.source_id=c2.source_id JOIN publications b ON a.publication_id=b.id JOIN spectral_types c ON a.source_id=c.source_id JOIN sources d ON a.source_id=d.id WHERE c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND a.instrument_id=6 AND (c.spectral_type<=15.5 AND c.spectral_type>=10.0) AND ((color>=1.33 AND (c.spectral_type=10 OR c.spectral_type=10.5)) OR (color>=1.33 AND (c.spectral_type=11 OR c.spectral_type=11.5)) OR (color>=1.67 AND c.spectral_type=12) OR (color>=1.63 AND (c.spectral_type=13 OR c.spectral_type=13.5)) OR (color>=1.86 AND (c.spectral_type=14 OR c.spectral_type=14.5)) OR (color>=1.46 AND (c.spectral_type=15 OR c.spectral_type=15.5)))").fetchall()  
out = open('Files/dbquery_red_sp.txt', 'wb')                                                                     
cPickle.dump(result, out)                     

#quit()
#print result
#discovery ref
result = db.query.execute("SELECT a.id, a.designation, b.id, b.bibtex, b.shortname, c.gravity, c.spectral_type, (c1.magnitude-c2.magnitude) AS color FROM sources a JOIN photometry c1 ON a.id=c1.source_id JOIN photometry c2 ON a.id=c2.source_id JOIN publications b ON a.publication_id=b.id JOIN spectral_types c ON a.id=c.source_id WHERE c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND (c.spectral_type<=15.5 AND c.spectral_type>=10.0) AND ((color>=1.33 AND (c.spectral_type=10 OR c.spectral_type>=10.5)) OR (color>=1.33 AND (c.spectral_type=11 OR c.spectral_type=11.5)) OR (color>=1.67 AND (c.spectral_type=12 OR c.spectral_type=12.5)) OR (color>=1.63 AND (c.spectral_type=13 OR c.spectral_type=13.5)) OR (color>=1.86 AND (c.spectral_type=14 OR c.spectral_type=14)) OR (color>=1.46 AND (c.spectral_type=15 OR c.spectral_type=15.5)))").fetchall()  
out = open('Files/dbquery_red_dis.txt', 'wb')                                                                     
cPickle.dump(result, out)                     

quit()

# spex prism ref for standard
result = db.query.execute("SELECT a.source_id, a.filename, b.id, b.bibtex, b.shortname, c.gravity, c.spectral_type, (c1.magnitude-c2.magnitude) AS color FROM spectra a JOIN photometry c1 ON a.source_id=c1.source_id JOIN photometry c2 ON a.source_id=c2.source_id JOIN publications b ON a.publication_id=b.id JOIN spectral_types c ON a.source_id=c.source_id WHERE c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND a.instrument_id=6 AND (c.spectral_type<=15.5 AND c.spectral_type>=10.0) AND ((color>=1.32 AND color <=1.34 AND (c.spectral_type=10 OR c.spectral_type>=10.5)) OR (color>=1.32 AND color <=1.34 AND (c.spectral_type=11.5 OR c.spectral_type=11)) OR (color>=1.66 AND color<=1.68 AND c.spectral_type=12) OR (color>=1.62 AND color <=1.64 AND c.spectral_type=13) OR (color>=1.85 AND color <= 1.87 AND c.spectral_type=14) OR (color<=1.45 AND color>=1.47 AND (c.spectral_type=15 OR c.spectral_type=15.5)))").fetchall()  

out = open('Files/dbquery_st_sp.txt', 'wb')
cPickle.dump(result, out)
#asciitable.write(result, 'Files/dbquery_lg.txt')
#quit()

# discovery ref for spec standard
result = db.query.execute("SELECT a.designation, b.id, a.id, b.bibtex, b.shortname, c.gravity, c.spectral_type, (c1.magnitude-c2.magnitude) AS color FROM sources a JOIN photometry c1 ON a.id=c1.source_id JOIN photometry c2 ON a.id=c2.source_id JOIN publications b ON a.publication_id=b.id JOIN spectral_types c ON a.id=c.source_id WHERE (c.spectral_type<=15.5 AND c.spectral_type>=10.0) AND c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND ((color>=1.32 AND color <=1.34 AND (c.spectral_type>=10.5 OR c.spectral_type=10)) OR (color>=1.32 AND color <=1.34 AND (c.spectral_type=11 OR c.spectral_type=11.5)) OR (color>=1.66 AND color<=1.68 AND (c.spectral_type=12 OR c.spectral_type=12.5)) OR (color>=1.62 AND color <=1.64 AND (c.spectral_type=13 OR c.spectral_type=13.5)) OR (color>=1.85 AND color <= 1.87 AND (c.spectral_type=14 OR c.spectral_type=14.5)) OR (color<=1.45 AND color>=1.47 AND (c.spectral_type=15 OR c.spectral_type=15.5)))").fetchall()  

out = open('Files/dbquery_st.txt', 'wb')
cPickle.dump(result, out)
#asciitable.write(result, 'Files/dbquery_lg.txt')
quit()



result = db.query.execute("SELECT a.id, a.names, b.spectral_type, b.gravity, (c1.magnitude-c2.magnitude) AS color, d.wavelength, d.flux, d.unc FROM sources a JOIN spectral_types b ON a.id=b.source_id JOIN photometry c1 ON a.id=c1.source_id JOIN photometry c2 ON a.id=c2.source_id JOIN spectra d ON a.id=d.source_id WHERE b.gravity IS NULL AND c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND d.instrument_id=6 AND b.regime='OPT' AND ((color>=1.3 AND b.spectral_type=10) OR (color>=1.35 AND b.spectral_type=11) OR (color>=1.48 AND b.spectral_type=12) OR (color>=1.64 AND b.spectral_type=13) OR (color>=1.69 AND b.spectral_type=14) OR (color>=1.72 AND b.spectral_type=15))").fetchall()

asciitable.write(result, 'Files/dbquery_fg1.txt')
quit()

result = db.query.execute("SELECT DISTINCT a.id, a.names, b.spectral_type, b.gravity, (c1.magnitude-c2.magnitude) AS color FROM sources a JOIN spectral_types b ON a.id=b.source_id JOIN photometry c1 ON a.id=c1.source_id JOIN photometry c2 ON a.id=c2.source_id JOIN spectra d ON a.id=d.source_id WHERE b.gravity IS NULL AND c1.band='J' AND c1.system=2 AND c2.band='Ks' AND c2.system=2 AND (SELECT COUNT(*) FROM spectra d WHERE d.source_id=a.id AND d.instrument_id=6)=0 AND b.regime='OPT' AND ((color>=1.3 AND b.spectral_type=10) OR (color>=1.35 AND b.spectral_type=11) OR (color>=1.48 AND b.spectral_type=12) OR (color>=1.64 AND b.spectral_type=13) OR (color>=1.69 AND b.spectral_type=14) OR (color>=1.72 AND b.spectral_type=15))").fetchall()

#print result
asciitable.write(result, 'Files/dbquery_fg2.txt')
