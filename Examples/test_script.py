from src.Main import MESMER_API
import src.Visualisation as Vis


me = MESMER_API()
me.parse_me_xml('MFfit_out.xml')
#Vis.plotPES(me)
Vis.plot_species_profile(me,Temperature=429,Species=["OOCO[C]=O","COC=O", "O=C=O"])
print('done')