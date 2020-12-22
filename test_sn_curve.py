from sn_curve_iso_6336 import SN_curve_ISO_6336

def check_results():
    material = SN_curve_ISO_6336('18CrNiMo6', 1e3, 3e6,2520,1050,1e5,5e7,2400,1550)

    assert material.name == '18CrNiMo6'
    p_F, p_H = material.calc_slope()[:2]
    assert 9.1 < p_F < 9.2
    assert 14.1 < p_H < 14.3
    for name, val in material.as_dict().items():
        print(name, val)
    print('Great job, you made it to the end of the code checks!')

    print(repr(material))
    material.plot_SN_curve()

check_results()