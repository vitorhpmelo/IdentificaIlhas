#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class bar():
    def __init__(self, id, type, contador):
        self.id = id
        self.type = type
        self.i = contador
        self.V = 1
        self.teta = 0
        self.Sbase = 100
        self.Vbase = 138
        self.Pg = 0
        self.Qg = 0
        self.Pd = 0
        self.Qd = 0
        self.Bs = 0
        self.nloads = 0
        self.nshunts = 0
        self.ngds = 0


class branch():
    def __init__(self, id, de, para, type, i):
        self.id = id
        self.de = de
        self.para = para
        self.type = type
        self.i = i
        self.x = -1
        self.r = -1
        self.ykm = 0
        self.Y = np.zeros((2, 2), dtype=complex)
        self.bsh = -1
        self.tap = -1
        self.limPA = -999
        self.flagLimP = 0

    def cykm(self):
        self.ykm = 1/complex(self.r, self.x)

    def twoPortCircuit(self):
        if self.type == 1:
            self.Y[0][0] = self.ykm+self.bsh
            self.Y[1][1] = self.ykm+self.bsh
            self.Y[1][0] = -self.ykm
            self.Y[0][1] = -self.ykm
        elif self.type == 2:
            self.Y[0][0] = ((1/self.tap)**2)*self.ykm
            self.Y[1][1] = self.ykm
            self.Y[1][0] = -(1/self.tap)*self.ykm
            self.Y[0][1] = -(1/self.tap)*self.ykm

    def Pf(self, grafo, flagT):
        k = self.de
        m = self.para
        if flagT == 0:
            P = (grafo[k].V**2)*np.real(self.Y[0][0]) + grafo[k].V * grafo[m].V*(
                np.real(self.Y[0][1])*np.cos(grafo[k].teta-grafo[m].teta)
                + np.imag(self.Y[0][1])*np.sin(grafo[k].teta-grafo[m].teta))
            return P
        elif flagT == 1:
            P = (grafo[m].V**2)*np.real(self.Y[1][1]) + grafo[m].V * grafo[k].V*(
                np.real(self.Y[1][0])*np.cos(grafo[m].teta-grafo[k].teta)
                + np.imag(self.Y[0][1])*np.sin(grafo[m].teta-grafo[k].teta))
            return P
        else:
            return 0

    def Qf(self, grafo, flagT):
        k = self.de
        m = self.para
        if flagT == 0:
            Qf = -(grafo[k].V**2)*np.imag(self.Y[0][0]) - grafo[k].V * grafo[m].V*(
                np.imag(self.Y[0][1])*np.cos(grafo[k].teta-grafo[m].teta)
                - np.real(self.Y[0][1])*np.sin(grafo[k].teta-grafo[m].teta))
            return Qf
        elif flagT == 1:
            Qf = -(grafo[m].V**2)*np.imag(self.Y[1][1]) - grafo[m].V * grafo[k].V*(
                np.imag(self.Y[1][0])*np.cos(grafo[m].teta-grafo[k].teta)
                - np.real(self.Y[1][0])*np.sin(grafo[m].teta-grafo[k].teta))
            return Qf
        else:
            return 0

    def dPfdt(self, grafo, flagT, var):
        k = self.de
        m = self.para
        Vk = grafo[k].V
        Vm = grafo[m].V
        tk = grafo[k].teta
        tm = grafo[m].teta
        if flagT == 0:
            Bkm = np.imag(self.Y[0][1])
            Gkm = np.real(self.Y[0][1])
            if k == var:  # dPkm/dtk
                return Vk*Vm*(Bkm*np.cos(tk-tm)-Gkm*np.sin(tk-tm))
            elif m == var:  # dPkm/dtm
                return Vk*Vm*(-Bkm*np.cos(tk-tm)+Gkm*np.sin(tk-tm))
            else:
                return 0
        elif flagT == 1:
            Bmk = np.imag(self.Y[1][0])
            Gmk = np.real(self.Y[1][0])
            if k == var:  # dPmk/dtk
                return Vk*Vm*(-Bmk*np.cos(tk-tm)-Gmk*np.sin(tk-tm))
            elif m == var:  # dPmk/dtm
                return Vk*Vm*(Bmk*np.cos(tk-tm)+Gmk*np.sin(tk-tm))
            else:
                return 0

    def dPfdV(self, grafo, flagT, var):
        k = self.de
        m = self.para
        Vk = grafo[k].V
        Vm = grafo[m].V
        tk = grafo[k].teta
        tm = grafo[m].teta
        if flagT == 0:  # dPkm
            Gkk = np.real(self.Y[0][0])
            Bkm = np.imag(self.Y[0][1])
            Gkm = np.real(self.Y[0][1])
            if k == var:  # dPkm/dVk
                return 2*Gkk*Vk + Vm*(Bkm*np.sin(tk-tm)+Gkm*np.cos(tk-tm))
            elif m == var:
                return Vk*(Bkm*np.sin(tk-tm)+Gkm*np.cos(tk-tm))
            else:
                return 0
        elif flagT == 1:
            Gmm = np.real(self.Y[1][1])
            Bmk = np.imag(self.Y[1][0])
            Gmk = np.real(self.Y[1][0])
            if k == var:  # dPkm/dVk
                return Vm*(-Bmk*np.sin(tk-tm)+Gmk*np.cos(tk-tm))
            elif m == var:  # dPkm/dVm
                return 2*Gmm*Vm + Vk*(-Bmk*np.sin(tk-tm)+Gmk*np.cos(tk-tm))
            else:
                return 0

    def dQfdt(self, grafo, flagT, var):
        k = self.de
        m = self.para
        Vk = grafo[k].V
        Vm = grafo[m].V
        tk = grafo[k].teta
        tm = grafo[m].teta
        if flagT == 0:  # dQkm
            Bkm = np.imag(self.Y[0][1])
            Gkm = np.real(self.Y[0][1])
            if k == var:  # dQkm/dtk
                return Vk*Vm*(Bkm*np.sin(tk-tm)+Gkm*np.cos(tk-tm))
            elif m == var:  # dQkm/dtm
                return Vk*Vm*(-Bkm*np.sin(tk-tm)-Gkm*np.cos(tk-tm))
        elif flagT == 1:  # dQmk
            Bmk = np.imag(self.Y[1][0])
            Gmk = np.real(self.Y[1][0])
            if k == var:  # dQmk/dtk
                return Vk*Vm*(Bmk*np.sin(tk-tm)-Gmk*np.cos(tk-tm))
            elif m == var:  # dQmk/dtm
                return Vk*Vm*(-Bmk*np.sin(tk-tm)+Gmk*np.cos(tk-tm))
            else:
                return 0

    def dQfdV(self, grafo, flagT, var):
        k = self.de
        m = self.para
        Vk = grafo[k].V
        Vm = grafo[m].V
        tk = grafo[k].teta
        tm = grafo[m].teta
        if flagT == 0:  # dQkm
            Bkk = np.imag(self.Y[0][0])
            Bkm = np.imag(self.Y[0][1])
            Gkm = np.real(self.Y[0][1])
            if k == var:  # dQkm/dvk
                return -2*Bkk*Vk+Vm*(-Bkm*np.cos(tk-tm)+Gkm*np.sin(tk-tm))
            elif m == var:  # dQkm/dvm
                return Vk*(-Bkm*np.cos(tk-tm)+Gkm*np.sin(tk-tm))
            else:
                return 0
        elif flagT == 1:
            Bmm = np.imag(self.Y[1][1])
            Bmk = np.imag(self.Y[1][0])
            Gmk = np.real(self.Y[1][0])
            if k == var:  # dQmk/dvk
                return Vm*(-Bmk*np.cos(tk-tm)-Gmk*np.sin(tk-tm))
            elif m == var:  # dQmk/dvm
                return -2*Bmm*Vm + Vk*(-Bmk*np.cos(tk-tm)-Gmk*np.sin(tk-tm))
            else:
                return 0


class node_graph():
    def __init__(self, id, bar):
        self.V = 1
        self.teta = 0
        self.Bs = 0
        self.adjk = dict()
        self.adjm = dict()
        self.ladjk = []
        self.ladjm = []
        self.id = id
        self.bar = bar
        self.V = 1
        self.teta = 0
        self.FlagBS = 0
        self.Bs = 0
        self.adjk = dict()
        self.adjm = dict()
        self.xy=None

    def P(self, graph):
        P = 0
        for key, item in self.adjk.items():
            P = P+item.Pf(graph, 0)
        for key, item in self.adjm.items():
            P = P+item.Pf(graph, 1)
        if np.abs(P) < 1e-12:
            P = 0
        return P

    def Q(self, graph):
        if self.FlagBS == 0:
            Q = 0
        else:
            Q = -self.Bs*self.V**2
        for key, item in self.adjk.items():
            Q = Q+item.Qf(graph, 0)
        for key, item in self.adjm.items():
            Q = Q+item.Qf(graph, 1)
        if np.abs(Q) < 1e-12:
            Q = 0
        return Q

    def dPdt(self, graph, bar):
        if self.i == bar:
            dPdt = 0
            for key, item in self.adjk.items():
                dPdt = dPdt+item.dPfdt(graph, FlagT=0, var=bar)
            for key, item in self.adjm.items():
                dPdt = dPdt+item.dPfdt(graph, FlagT=0, var=bar)
            return dPdt
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdt(graph, 0, bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdt(graph, 1, bar)
        else:
            return 0

    def dPdV(self, graph, bar):
        dPdV = 0
        if self.i == bar:
            for key, item in self.adjk.items():
                dPdV = dPdV+item.dPfdV(graph, FlagT=0, var=bar)
            for key, item in self.adjm.items():
                dPdV = dPdV+item.dPfdV(graph, FlagT=0, var=bar)
            return dPdV
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdV(graph, 0, bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dPfdV(graph, 1, bar)
        else:
            return 0

    def dQdt(self, graph, bar):
        dQdt = 0
        if self.i == bar:
            for key, item in self.adjk.items():
                dQdt = dQdt+item.dQdt(graph, FlagT=0, var=bar)
            for key, item in self.adjm.items():
                dQdt = dQdt+item.dQdt(graph, FlagT=0, var=bar)
            return dQdt
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dQdt(graph, 0, bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dQdt(graph, 1, bar)
        else:
            return 0

    def dQdV(self, graph, bar):
        if self.i == bar:
            if self.FlagBS == 0:
                dQdV = 0
            else:
                dQdV = -2*self.Bs*self.V
            for key, item in self.adjk.items():
                dQdV = dQdV+item.dQdV(graph, FlagT=0, var=bar)
            for key, item in self.adjm.items():
                dQdV = dQdV+item.dQdV(graph, FlagT=0, var=bar)
            return dQdV
        elif str(self.i)+"-"+str(bar) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dQdV(graph, 0, bar)
        elif str(bar)+"-"+str(self.i) in self.adjk.keys():
            return self.adjk[str(self.i)+"-"+str(bar)].dQdV(graph, 1, bar)
        else:
            return 0


class medida():
    """
    objeto do tipo medida que constitui o plano de medições, tem como atributos
    @param: instalado 1 medida instalada, 2 PMU instalada, 0 Sem PMU ()
    @param: tipo 0 injeção, 2 fluxo, 5 MFS_V, 6 MFS_I
    @param: id barra de
    @param: id barra para
    """

    def __init__(self, instalado, tipo, de, para):
        self.instalado = instalado
        self.tipo = tipo
        self.de = de
        self.para = para
