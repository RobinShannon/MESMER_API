from MESMER_API.src.Main import MESMER_API
import MESMER_API.src.Visualisation as Vis


me = MESMER_API()
me.parse_me_xml('mestemplate.xml')
Vis.plotPES(me)
print('done')