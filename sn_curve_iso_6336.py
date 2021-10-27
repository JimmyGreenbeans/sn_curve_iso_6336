import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (ScalarFormatter, MultipleLocator)
import math
import pandas as pd
import time
from typing import Tuple


class SnCurveIso6336:
    '''
    Represents individual material data for gear calculations according to 
    ISO 6336
    
    '''
    
    def __init__(self, name: str, N_F_stat: int, N_F_d: int, sig_FP_stat: float, sig_FE: float, 
                N_H_stat: int, N_H_d: int, sig_HP_stat: float, sig_H_lim: float, 
                lim_pit_perm: bool = False, red_life_fac: bool = False):
        """
        Args:
            name (str): material name
            N_F_stat (int): number of load cycles for endurance limit of tooth root
            N_F_d (int): upper number of load cycles for static limit of tooth root
            sig_FP_stat (float): static limit of tooth root
            sig_FE (float): endurance limit of tooth root. notch sensitivity factor of
                            the standard reference test gear must be included, 
                            sig_FE = sig_F_lim * Y_ST
            N_H_stat (int): number of load cycles for endurance limit of tooth flank
            N_H_d (int): upper number of load cycles for static limit of tooth flank
            sig_HP_stat (float): static limit of tooth flank
            sig_H_lim (float): endurance limit of tooth flank
            lim_pit_perm (bool, optional): Is limited pitting according to ISO 6336-2 
                                            permitted? Defaults to False.
            red_life_fac (bool, optional): Reduction of life factors Y_NT and Z_NT to 
                                            0,85 according to ISO 6336 from N_F_d and 
                                            N_H_d to 10^10 load cycles? Defaults to 
                                            False.
        """

        self.name = name
        self.N_F_stat = N_F_stat
        self.N_F_d = N_F_d
        self.sig_FP_stat = sig_FP_stat
        self.sig_FE = sig_FE
        self.N_H_stat = N_H_stat
        self.N_H_d = N_H_d
        self.sig_HP_stat = sig_HP_stat
        self.sig_H_lim = sig_H_lim
        self.lim_pit_perm = lim_pit_perm
        self.red_life_fac = red_life_fac
        self._my_var = 'akgflkaga'
        self._radius = 10.

        # MS: - Sind diese Attribute wirklich alle public? 
        #       Wenn ja, besser als Property mit validierung, oder zumindest Dockstrings inkl type für jedes Attribut schreiben
        #       Wenn nein, als private self._var_name benennen
        #     - init-parameter auf Type und Intervall prüfen, bei Fehlschlägen Exception auslösen
        #       Gilt ebenfalls für alle anderen Public-Methoden
        #     - Welche Methoden müssen Public sein?
        #     - __str__ und __repr__ Methoden implementieren
        

    # region properties
    
    @property
    def my_var(self) -> str: 
        '''Liefert die variable tralala oder legt sie fest'''
        return self._my_var
    @my_var.setter
    def my_var(self, value):
        self._my_var = value
        print("my_var gesetzt!")

    @property
    def CircleArea(self) -> float: 
        '''Liefert die Flaeche des Kreises mit Radius R'''
        return self._radius**2 * 3.14
    # endregion

    def __repr__(self):
        return f'SnCurveIso6336({self.name},{self.N_F_stat},{self.N_F_d})'

    def calc_slope(self) -> Tuple[float]:
        """
        Calculates the exponent (slope) of the S-N curve for the tooth root and the flank

        Returned value p_H for pitting is for stress; to convert for torque, these value must be divided by 2

        Returns:
            tuple[float]: p_F, p_H, p_H_lim_pit, p_F_red_life_fac, p_H_red_life_fac
        """
        
        # MS: Methoden ohne Parameter und "billigen berechnungen" besser als ReadOnly Property und Name 
        #     ohne Imperativ, also nur slope
        # Tipp: - abkürzung log = np.log und np. in formeln weg lassen
        #       - self-attribute in loc variablen schreiben und self. in formeln weg lassen
        #       - Leerzeichen vor und nach mathematischen operatoren. immer!  
        p_F = (np.log(self.N_F_d)-np.log(self.N_F_stat))/(np.log(self.sig_FP_stat)-np.log(self.sig_FE))
        p_H = (np.log(self.N_H_d)-np.log(self.N_H_stat))/(np.log(self.sig_HP_stat)-np.log(self.sig_H_lim))

        # bsp mit Abkürzungen
        # p_F = (log(N_F_d) - log(N_F_stat)) / (log(sig_FP_stat) - log(sig_FE))
        # p_H = (log(N_H_d) - log(N_H_stat)) / (log(sig_HP_stat) - log(sig_H_lim))
        
        # MS: == 1 -> Kein Vergleich mit ints. lim_pit_perm ist bereits ein bool. Besser: if self.lim_pit_perm: '
        if self.lim_pit_perm == 1: # when limited pitting is permitted s-n curve for flank has two different slopes
            p_H = (np.log(1e7)-np.log(self.N_H_stat))/(np.log(self.sig_HP_stat)-np.log(0.5*(self.sig_H_lim+self.sig_HP_stat)))
            p_H_lim_pit = (np.log(1e9)-np.log(1e7))/(np.log(0.5*(self.sig_H_lim+self.sig_HP_stat))-np.log(self.sig_H_lim))
        else:
            p_H_lim_pit = '-' # MS: warum ein string? besser: -1. Kann gut geprüft werden, ist schneller und das Rückgabetupel bleibt numerisch
        
        if self.red_life_fac == 1:
            p_F_red_life_fac = (np.log(1e10)-np.log(self.N_F_d))/(np.log(self.sig_FE)-np.log(0.85*self.sig_FE))
            p_H_red_life_fac = (np.log(1e10)-np.log(self.N_H_d))/(np.log(self.sig_H_lim)-np.log(0.85*self.sig_H_lim))
        else:
            p_F_red_life_fac = '-' # MS: warum ein string? besser: -1. Kann gut geprüft werden, ist schneller und das Rückgabetupel bleibt numerisch
            p_H_red_life_fac = '-' # MS: warum ein string? besser: -1. Kann gut geprüft werden, ist schneller und das Rückgabetupel bleibt numerisch
                   
        return p_F, p_H, p_H_lim_pit, p_F_red_life_fac, p_H_red_life_fac
    
    
    def as_dict(self) -> dict:
        # MS: Convention bei solchen methoden: to_dict()
        """
        Transforms material entry to a dictionary, used in function "export_data"'''

        Returns:
            dict: dictionary with following keys:
                  'name', 'N_F_stat', 'N_F_d', 'sig_FP_stat', 'sig_FE', 'N_H_stat', 
                  'N_H_d', 'sig_HP_stat', 'sig_H_lim', 'lim_pit_perm', 'red_life_fac', 
                  'p_F', 'p_H', 'p_H_lim_pit', 'p_F_red_life_fac', 'p_H_red_life_fac'
        """

        p_F, p_H, p_H_lim_pit, p_F_red_life_fac, p_H_red_life_fac = self.calc_slope()

        return {'name': self.name, 'N_F_stat': self.N_F_stat, 'N_F_d': self.N_F_d, 
                'sig_FP_stat': self.sig_FP_stat, 'sig_FE': self.sig_FE, 
                'N_H_stat': self.N_H_stat, 'N_H_d': self.N_H_d, 
                'sig_HP_stat': self.sig_HP_stat, 'sig_H_lim': self.sig_H_lim, 
                'lim_pit_perm': self.lim_pit_perm, 'red_life_fac': self.red_life_fac, 
                'p_F': p_F, 'p_H': p_H, 'p_H_lim_pit': p_H_lim_pit, 
                'p_F_red_life_fac': p_F_red_life_fac, 'p_H_red_life_fac':p_H_red_life_fac}
    
    
    def plot_SN_curve(self, save_file: bool = False):
        '''
        plots the S-N curves for the tooth root and the flank of a S-N curve object    
        not included in permissable stress are the factors according to ISO 6336 like surface factor
        '''

        # speichern als PDF soll optional gehen
        self.sig_FE_red_life_fac = self.sig_FE if self.red_life_fac == 0 else 0.85 * self.sig_FE
        
        N = [1, self.N_F_stat, self.N_F_d, 1e10]
        stress = [self.sig_FP_stat, self.sig_FP_stat, self.sig_FE, self.sig_FE_red_life_fac]
        plt.loglog(N, stress, '-')
        
        N_lim_perm = self.N_H_stat if self.lim_pit_perm == 0 else 1e7
        sig_H_lim_perm = self.sig_HP_stat if self.lim_pit_perm == 0 else 0.5*(self.sig_H_lim+self.sig_HP_stat)
        self.sig_H_lim_red_life_fac = self.sig_H_lim if self.red_life_fac ==0 else 0.85 * self.sig_H_lim
        
        N = [1, self.N_H_stat, N_lim_perm, self.N_H_d, 1e10]
        stress = [self.sig_HP_stat, self.sig_HP_stat, sig_H_lim_perm, self.sig_H_lim, self.sig_H_lim_red_life_fac]

        plt.loglog(N, stress, '-')
    
        all_values = list(set([self.sig_FP_stat, self.sig_FE, self.sig_HP_stat, self.sig_H_lim, self.sig_FE_red_life_fac]))
        min_value = math.floor(min(all_values)/100.)*100
        max_value = math.ceil(max(all_values)/100.)*100
        
        plt.xlabel('Number of load cycles')
        plt.ylabel('Permissable bending stress in MPa')
        plt.title('S-N curves for foot and flank of ' + self.name)
        plt.legend(['foot', 'flank'])
        
        # get the current axes, creating them if necessary:
        ax = plt.gca()
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(ScalarFormatter())
        ax.set_ylim(min_value-100, max_value+100)
        ax.yaxis.set_major_locator(MultipleLocator(200))
        ax.tick_params(axis='y', which='minor', labelbottom=False, labeltop=False, labelleft=False, labelright=False) # Bug in Matplotlib, sonst Werte auf y-Achse doppelt
        ax.grid(True, which='major',linewidth=0.5)
        ax.grid(True, which='minor', linestyle='--', linewidth=0.3)
        ax.set_xlim(self.N_F_stat/10, 1e10)
        
        ax2 = ax.twinx()
        ax2.set_yscale('log')
        ax2.set_ylim(min_value-100, max_value+100)
        ax2.yaxis.set_major_formatter(ScalarFormatter())
        ax2.yaxis.set_minor_formatter(ScalarFormatter())
        ax2.set_yticks(all_values)
        ax2.tick_params(axis='y', which='minor', labelbottom=False, labeltop=False, labelleft=False, labelright=False) # Bug in Matplotlib, sonst Werte auf y-Achse doppelt
        ax2.grid(True, which='major', linestyle='--', linewidth=0.3)
        
        #plt.xlim(1e2, 1e10)
        if save_file:
            plt.savefig('sn_curve_%s.pdf' % (self.name), format='pdf') 
        plt.show()

    
    def calc_sig_perm(self, load_cycles: int) -> Tuple[float]:
        """
        returns permissible stress for foot and flank for given number of load cycles

        Args:
            load_cycles (int): Number of stress cycles for which the perm. stresses should be returned

        Returns:
            Tuple[float]: perm stress foot, perm stress flank
        """
        
        p_F, p_H, p_H_lim_pit, p_F_red_life_fac, p_H_red_life_fac = self.calc_slope()
        
        # Calculation of perm. stress sig_perm_F
        if load_cycles <= self.N_F_stat: 
            sig_perm_F = self.sig_FP_stat
        elif load_cycles >= self.N_F_d:
            if self.red_life_fac:
                if load_cycles > 1e10:
                    sig_perm_F = self.sig_FE * 0.85
                else:
                    sig_perm_F = self.sig_FE * 0.85 * (1e10/load_cycles)**(1/p_F_red_life_fac)
            else:
                sig_perm_F = self.sig_FE
        elif load_cycles > self.N_F_stat and load_cycles < self.N_F_d:
            sig_perm_F = self.sig_FE * (self.N_F_d/load_cycles)**(1/p_F)
        
        # Calculation of perm. stress sig_perm_H
        if load_cycles <= self.N_H_stat: 
            sig_perm_H = self.sig_HP_stat
        elif load_cycles > self.N_H_stat and load_cycles < self.N_H_d:
            if self.lim_pit_perm:
                if load_cycles <= 1e7:
                    sig_perm_H = (self.sig_H_lim + self.sig_HP_stat) / 2 * (1e7/load_cycles)**(1/p_H)
                else:
                    sig_perm_H = self.sig_H_lim * (self.N_H_d/load_cycles)**(1/p_H_lim_pit)
            else:
                sig_perm_H = self.sig_H_lim * (self.N_H_d/load_cycles)**(1/p_H)
        elif load_cycles >= self.N_H_d:
            if self.red_life_fac:
                if load_cycles > 1e10:
                    sig_perm_H = self.sig_H_lim * 0.85
                else:
                    sig_perm_H = self.sig_H_lim * 0.85 * (1e10/load_cycles)**(1/p_H_red_life_fac)
            else:
                sig_perm_H = self.sig_H_lim
            
        #print('sig_perm_F: ' + str(sig_perm_F) + ', sig_perm_H: ' + str(sig_perm_H))
        return sig_perm_F, sig_perm_H
        
    
    def calc_N_perm(self, stress):
        '''returns permissible number of load cycles for given number of load cycles
        
        Args:
        stress (float): number of load cycles for which the permissible stress shall be calculated
        '''
        raise NotImplementedError("Sorry, time has running out.")
        # TODO
        #N_perm_F = self.N_F_d * (self.sig_FE/self.sig_FP_stat)**p
    
    
    def write_dat_file(self):
        '''
        creates and saves dat-file
        creates some interpolation points
        '''
        
        # create list of numbers from 100 to 10^99
        a1 = np.geomspace(1000, 1e10, 8)
        a2 = np.arange(1,10,1)
        NL = np.outer(a1, a2).flatten()
        NL = NL.tolist()
        NL.append(1e99)
        NL_mod = NL
        NL_mod.insert(0, 0e0)
        sig_perm_F = []
        sig_perm_H = []
        
        for number in NL_mod:
            foot, flank = self.calc_sig_perm(number)
            sig_perm_F.append(foot)
            sig_perm_H.append(flank)
        
        #inhalt = [sheet['A20':'D47'],sheet['F20':'I47'],sheet['K20':'N47'],sheet['P20':'S47'],sheet['U20':'X47'],sheet['Z20':'AC47'],sheet['AE20':'AH47'],sheet['AJ20':'AM47'],sheet['AO20':'AR47'],sheet['AT20':'AW47'],sheet['AY20':'BB47'],sheet['BD20':'BG47'],sheet['BI20':'BL47']]
        #inhalt = [list(i) for i in inhalt]
        #werkstoff = [sheet['B3'].value,sheet['G3'].value, sheet['L3'].value, sheet['Q3'].value, sheet['V3'].value, sheet['AA3'].value, sheet['AF3'].value, sheet['AK3'].value, sheet['AP3'].value, sheet['AU3'].value, sheet['AZ3'].value, sheet['BE3'].value, sheet['BJ3'].value] 
        
        # Schreiben der neuen Dateien   
        # Dateipfad und Name und Deklarierung als 'f'
        #with open("I:\\Technische_Berechnung\\Projekte\\AT\\Grundlagenuntersuchungen\\150427 Eigene Wöhlerlinien in KISSsoft\\Textdateien_Woehlerlinien_DNV_KISSsoft\\WL_" + "0" + "_" + self.name + ".dat", "w") as f:
        with open("WL_" + self.name + ".dat", "w") as f:
            f.write("-- -----------------------------------------------------------\n")
            f.write("-- File = WL_" + self.name + ".dat" + "\n")
            f.write("-- Erstellt am " + (time.strftime("%d/%m/%Y")) + " von Jürgen Hammele\n\n")
            f.write("-- Werte nach DNV Report ER-DE-ISO6336-04848-1 vom 02.05.2019\n")
            #f.write("-- Werkstoff " + werkstoff[i] + "\n")
            f.write("-- -----------------------------------------------------------\n\n")
            
            # Ausgabe von Schwingspielanzahl
            f.write("-- Data for significant no. of cycles (edge points for interpolation of Woehler line)\n")
            f.write(":TABLE FUNCTION EdgeCycle\n")
            f.write("\tINPUT X number TREAT NEXT_BIGGER\n")
            f.write("DATA\n\t")
            for i in range(len(NL)):
                f.write(str(i+1) + "\t")
            f.write("\n\t")
            for number in NL:
                f.write(str('{:.0e}'.format(number)) + "\t")
            f.write("\nEND\n\n")
            
            # Ausgabe von ertragbarer Spannung Flanke
            f.write("-- Data for Hertzian pressure sigH\n")
            f.write(":TABLE FUNCTION FlankSigH\n")
            f.write("\tINPUT X Cycles TREAT LOG\n")
            f.write("DATA\n\t")
            for number in NL_mod:
                f.write(str('{:.0e}'.format(number)) + "\t")
            f.write("\n\t")
            for sig in sig_perm_H:
                f.write(str('{:.1f}\t'.format(sig)))
            f.write("\nEND\n\n")
            
            # Ausgabe von ertragbarer Spannung Fuss
            f.write("-- Data for fatigue strength tooth root sigF\n")
            f.write(":TABLE FUNCTION FootSigF\n")
            f.write("\tINPUT X Cycles TREAT LOG\n")
            f.write("DATA\n\t")
            for number in NL_mod:
                f.write(str('{:.0e}'.format(number)) + "\t")
            f.write("\n\t")
            for sig in sig_perm_F:
                f.write(str('{:.1f}\t'.format(sig/2)))
            f.write("\nEND\n\n")


# create list with standard materials of COB
materials_COB = []
materials_COB.append(SnCurveIso6336('AT-01_18CrNiMo7-6(COB)_LN_190-3', 1e3,3e6,2520,1050,1e5,5e7,2400,1550,0,0))
materials_COB.append(SnCurveIso6336('AT-02_18CrNiMo7-6(COB)_LN_190-3_lim_pit_perm', 1e3,3e6,2520,1050,6e5,1e9,2400,1550,1,0))
materials_COB.append(SnCurveIso6336('AT-03_18CrNiMo7-6(COB)_LN_190-3_shot_peened_LN_523-1', 1e3,3e6,2520,1400,1e5,5e7,2400,1800,0,0))
materials_COB.append(SnCurveIso6336('AT-04_18CrNiMo7-6(COB)_LN_190-3_lim_pit_perm_shot_peened_LN_523-1', 1e3,3e6,2520,1400,6e5,1e9,2400,1800,1,0))
materials_COB.append(SnCurveIso6336('AT-05_20MnCr5(COB)_LN_190-2', 1e3,3e6,2520,1050,1e5,5e7,2400,1550,0,0))
materials_COB.append(SnCurveIso6336('AT-06_20MnCr5(COB)_LN_190-2_lim_pit_perm', 1e3,3e6,2520,1050,6e5,1e9,2400,1550,1,0))
materials_COB.append(SnCurveIso6336('AT-07_20MnCr5(COB)_LN_190-2_shot_peened_LN_523-1', 1e3,3e6,2520,1400,1e5,5e7,2400,1800,0,0))
materials_COB.append(SnCurveIso6336('AT-08_20MnCr5(COB)_LN_190-2_lim_pit_perm_shot_peened_LN_523-1', 1e3,3e6,2520,1400,6e5,1e9,2400,1800,1,0))
materials_COB.append(SnCurveIso6336('AT-09_42CrMo4(COB)_LN_194-1', 1e3,3e6,1440,900,1e5,2e6,1430,1100,0,0))
materials_COB.append(SnCurveIso6336('AT-10_18CrNiMo7-6(COB)_LN_190-3_shot_peened_LN_523-2', 1e3,3e6,2520,1250,1e5,5e7,2400,1700,0,0))
materials_COB.append(SnCurveIso6336('AT-11_18CrNiMo7-6(COB)_LN_190-3_lim_pit_perm_shot_peened_LN_523-2', 1e3,3e6,2520,1250,6e5,1e9,2400,1700,1,0))
materials_COB.append(SnCurveIso6336('AT-12_20MnCr5(COB)_LN_190-2_shot_peened_LN_523-2', 1e3,3e6,2520,1250,1e5,5e7,2400,1700,0,0))
materials_COB.append(SnCurveIso6336('AT-13_20MnCr5(COB)_LN_190-2 _lim_pit_perm_shot_peened_LN_523-2', 1e3,3e6,2520,1250,6e5,1e9,2400,1700,1,0))
materials_COB.append(SnCurveIso6336('AT-14_42CrMo4(COB)_LN 191-1_root_49_HRC_root_flank_56HRC', 1e3,3e6,1800,720,1e5,5e7,1952,1220,0,0))
materials_COB.append(SnCurveIso6336('AT-15_42CrMo4(COB)_LN 191-1_root_49_HRC_root_flank_56HRC_lim_pit_perm', 1e3,3e6,1800,720,6e5,1e9,1952,1220,1,0))


materials_ISO = []


def export_data(materials: list):
    """
    Exports a list of materials to an Excel file

    Args:
        materials (list): list of SN_curve_ISO_6336 objects
    """
    
    df = pd.DataFrame([x.as_dict() for x in materials])
    df.to_excel('out.xlsx', index=False)


def plot_SN_curve_flank(materials: list):
    '''plots the S-N curves in one plot for flank of all materials in a list of SN_curve_ISO_6336 objects
    
    not included in permissable stress are the factors according to ISO 6336 like surface factor
    
    Args:
        materials (list): list of SN_curve_ISO_6336 objects
        '''
    
    names = []
    S_stat_values = []
    S_D_values = []
    S_D_red_life_values = []
    
    for material in materials:
        N_lim_perm = material.N_H_stat if material.lim_pit_perm == 0 else 1e7
        sig_H_lim_perm = material.sig_HP_stat if material.lim_pit_perm == 0 else 0.5*(material.sig_H_lim+material.sig_HP_stat)
        material.sig_H_lim_red_life_fac = material.sig_H_lim if material.red_life_fac ==0 else 0.85 * material.sig_H_lim
        
        N = [1, material.N_H_stat, N_lim_perm, material.N_H_d, 1e10]
        stress = [material.sig_HP_stat, material.sig_HP_stat, sig_H_lim_perm, material.sig_H_lim, material.sig_H_lim_red_life_fac]

        plt.loglog(N, stress, '-')

        names.append(material.name)
        S_stat_values.append(material.sig_HP_stat)
        S_D_values.append(material.sig_H_lim)
        S_D_red_life_values.append(material.sig_H_lim_red_life_fac)
    all_values = list(set(S_stat_values + S_D_values + S_D_red_life_values))
    min_value = math.floor(min(all_values)/100.)*100
    max_value = math.ceil(max(all_values)/100.)*100
        
    plt.xlabel('Number of load cycles')
    plt.ylabel('Permissable contact stress in MPa')
    plt.title('S-N curve flank')
    plt.xlim(1e4, 1e10)
    plt.legend(names, prop={'size': 3})
    
    # get the current axes, creating them if necessary:
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    ax.set_ylim(min_value-100, max_value+100)
    ax.tick_params(axis='y', which='minor', labelbottom=False, labeltop=False, labelleft=False, labelright=False) # Bug in Matplotlib, sonst Werte auf y-Achse doppelt
    ax.grid(True, which='major',linewidth=0.5)
    ax.grid(True, which='minor', linestyle='--', linewidth=0.3)
    
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    ax2.set_ylim(min_value-100, max_value+100)
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    ax2.yaxis.set_minor_formatter(ScalarFormatter())
    ax2.set_yticks(all_values)
    ax2.tick_params(axis='y', which='minor', labelbottom=False, labeltop=False, labelleft=False, labelright=False) # Bug in Matplotlib, sonst Werte auf y-Achse doppelt
    ax2.grid(True, which='major', linestyle='--', linewidth=0.3)
        
    plt.savefig('sn_curve_flank.pdf', format='pdf')
    plt.show()
    

def plot_SN_curve_foot(materials: list):
    '''plots the S-N curves in one plot for foot of all materials in a list of SN_curve_ISO_6336 objects
    
    not included in permissable stress are the factors according to ISO 6336
    
    Args:
        materials (list): list of SN_curve_ISO_6336 objects
        '''
    
    names = []
    S_stat_values = []
    S_D_values = []
    S_D_red_life_values = []
    
    for material in materials:
        material.sig_FE_red_life_fac = material.sig_FE if material.red_life_fac ==0 else 0.85 * material.sig_FE
        
        N = [1, material.N_F_stat, material.N_F_d, 1e10]
        stress = [material.sig_FP_stat, material.sig_FP_stat, material.sig_FE, material.sig_FE_red_life_fac]

        plt.loglog(N, stress, '-')

        names.append(material.name)
        S_stat_values.append(material.sig_FP_stat)
        S_D_values.append(material.sig_FE)
        S_D_red_life_values.append(material.sig_FE_red_life_fac)
    all_values = list(set(S_stat_values + S_D_values + S_D_red_life_values))
    min_value = math.floor(min(all_values)/100.)*100
    max_value = math.ceil(max(all_values)/100.)*100
        
    plt.xlabel('Number of load cycles')
    plt.ylabel('Permissable bending stress in MPa')
    plt.title('S-N curve root')
    plt.xlim(1e3, 1e10)
    plt.legend(names, prop={'size': 3})
    
    # get the current axes, creating them if necessary:
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    ax.set_ylim(min_value-100, max_value+100)
    ax.yaxis.set_major_locator(MultipleLocator(200))
    ax.tick_params(axis='y', which='minor', labelbottom=False, labeltop=False, labelleft=False, labelright=False) # Bug in Matplotlib, sonst Werte auf y-Achse doppelt
    ax.grid(True, which='major',linewidth=0.5)
    ax.grid(True, which='minor', linestyle='--', linewidth=0.3)
    #ax.set_xlim(np.min(N_stat_werte)/10, 1e10)
    
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    ax2.set_ylim(min_value-100, max_value+100)
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    ax2.yaxis.set_minor_formatter(ScalarFormatter())
    ax2.set_yticks(all_values)
    ax2.tick_params(axis='y', which='minor', labelbottom=False, labeltop=False, labelleft=False, labelright=False) # Bug in Matplotlib, sonst Werte auf y-Achse doppelt
    ax2.grid(True, which='major', linestyle='--', linewidth=0.3)
        
    plt.savefig('sn_curve_foot.pdf', format='pdf')
    plt.show()

y = materials_COB[3]
y.write_dat_file()