classdef AllV
    properties (Constant)
        CurrentVersion = '1.16';
        TEMPLATESPACE=3;
        JAMESSPACE = 7;
        USERSPACE = 8;
end
    methods (Static)        

function version = Version(obj)
    version = obj.CurrentVersion;
end
[res, varargout] = AllVPcs(V, varargin)

                 result = PlotTriggers(Vall, varargin);
                  newT = CheckDoubleTrigger(DATA, id, varargin);
[dropi, gfit, details] = CalcTriggerDrop(DATA, clid,varargin)
    [score, details] = CalcFitScore(fit);
                  DATA = CompareSpikeLists(DATA, ida, idb, iida, iidb)
       [fit, details] = CombineFits(DATA, fits,varargin)
                  score = CheckEllipse(DATA, C, varargin)
                             CheckLogs(dir, varargin);
  [bestfit, details] = ChooseBestFit(fits, varargin)
                 DATA = SetExpt(DATA);
                            AddCompareMenu(DATA, varargin);
                            MakeProbeMenu(DATA, sm, varargin);
                DATA = ApplyRefCut(a,b, varargin);
                [bestspace, DATA] = FindBestSpace(DATA, C, varargin);
               value = GetValue(DATA, type, varargin);
                DATA = PlotBestSpace(DATA, varargin);
                DATA = FinishNewData(DATA, varargin);
[xyr, ssq] = InvertEllipse(E, x,y)
DATA = SetData(DATA, op, varargin);
args = NewArgs(DATA, type);
Cd = UpdateClusterDetails(DATA, varargin);
CheckOneXY(DATA,varargin);
C = CompareMeans(Ca, Cb, varargin);
need = NeedEV(C);

[id, DATA] = Retrigger(DATA, varargin)
           fit = SetFit(fit, type, varargin)

yes = NeedNDspace(C, varargin)
 probes = GetProbeList(DATA);
          dv = GetdVdY(DATA,varargin);
      DATA = SetTemplatePlots(DATA, varargin)
[DATA, C, fit] = SetClusterFromFit(DATA, C, fit, varargin)

 [p, pid] = GetProbeNumber(DATA, chspk, n, varargin)

 DATA = CheckTrigTimes(DATA,varargin)
 [needc, C]  = NeedToQuantify(C, fullvname, varargin);
 
 DATA = ReadConfig(DATA, name, varargin)
 
 FullVButtonPressed(a,b)
 
 FullVButtonReleased(a,b)
 ScrollFullV(a,b)
 DATA = mysetappdata(DATA, str, val)

 val = mygetappdata(DATA, str)
[Lratio,iso_distance] = CalcLRatio(X,gmm_fit,comp_idx,cluster_labels, varargin)

 [AllVoltages, DATA] = BuildAllV(DATA, id, spts,varargin)

 val = myisappdata(DATA,str)

 DATA = CloseLog(DATA)

 C = CheckClusterFields(C, varargin)

 [V, DATA] = ReadSpikeFiles(DATA, name)

 str = IDStr(DATA, varargin)

 C = CheckScoreScaling(DATA, C)
 
 DATA = CalcIsolation(DATA, cnum, varargin)

 ClusterFromPoints(DATA)

 type = WhichPlotType(C,clnum, plottype)

 [DataClusters, FullVData] = LoadDataClusters(DATA)

 [DATA, Vall, ispk, newdata] = SetupVall(DATA,Vall, ispk, newdata)

  [true, each] =  NeedTemplateForCluster(C, all)

 DATA = GetGuiState(DATA, F)

 DATA = ResetDataForNewProbe(DATA)

 DATA = ReadFromLog(DATA)

 sz = memsize(X)

 S = SmallCluster(C)
 str = SpaceLabels(DATA, C, space);
 
 [true, cluster] = ClusterIsSet(C, cl)

 [DATA, ClusterDetails] = LoadClusterDetails(DATA, varargin)

 [DATA, DataClusters, success] = LoadTrigTimes(DATA, checktimes, varargin)

 [id,  th, details] = TriggerV(DATA, rV)

 [res, C,D] = AutoCutOne(DATA, exname, ispk, args)

 res = AutoCutAll(ispk, toplevel, Vall,DATA, args)

 first = PrevCluster(C)

 res = QuantifyQuickClusters(DATA, ispk, varargin) 

 need = NeedMore(C, Evec)

 C = CheckForMean(DATA,C)

 DATA = ReClassify(DATA, varargin)

 [DATA, Template] = CalcTemplateScores(DATA,  varargin)

 [DATA, Template] = CalcTemplatesFromMean(DATA, MeanSpike, varargin);

 tt = TimeMark(tt, str, show)

 good = GoodCluster(C)

 it = findobj(DATA, name);
 
 C = ClusterFromBoundary(E, C)

 E = BoundaryFromCluster(E, C, n)

 DATA = NextPCs(DATA)

 DATA = SetPCs(DATA, replot, reapply)

 DATA = CalcICA(DATA,ns)    

 [C, Evec, pcs, dip, chspk, errs, details] = CalcPCs(DATA, AllVoltages, nprobepc,  varargin)

 OldPlotFullV(V, ispk, DATA)

 distance  = gmdistance(G)

 [theta, c, details] = BestGMAngle(x,y, test, varargin)

 [theta, c, details] = BestAngleGM(xy, G, dipfit, varargin)

 [theta, c, details] = BestAngle(x,y, test, varargin)

 DATA = CheckClusterMarks(Clusters, DATA)

 CheckClusters(Clusters, str, varargin)

 SaveMeanSpikeOnly(DATA, outname) 

 CheckClusterValues(DATA, C)

 [DATA, id] = SaveClusters(DATA, outname,varargin)

 C = StripClusters(Clusters)

 name = SpkFileName(DATA, varargin)

 probe = ProbeNumber(DATA)

 SaveSpikes(DATA, id, name)

 tcut = IsTemplateCut(E)    

 [E, cluster] = CutAndSave(DATA, varargin)

 [distance, obj, xy, details] = TemplateSpace(DATA, varargin)

 DATA = TemplateGMFits(DATA)        

 [distance, obj, xy, details] = BestSpace(DATA, varargin)

 [G, D, all] = GMfit(X, nd, nr, varargin)

 [G, xy] = IterateTemplateFit(DATA, G)

 [xy, cid, details] = ProjectND(DATA, best, obj)

 Spikes = MakeJamesSpikes(DATA)

 [DATA, details]  = JamesAutoCut(DATA, varargin)

 [DATA, details]  = EckerAutoCut(DATA, varargin)

 [E, Scores, tbs, xy, details]  = AutoCut(DATA, varargin)

 DATA = AddErr(DATA,varargin)

 C = CutAndPlot(x,y, energy)

 C = OptimizeVarE(DATA)

 [DATA, E] = OptimizeBoundary(DATA);

 [C, details] = OptimizeClusterBoundary(DATA)

 FullV = SetFullVNames(DATA, FullV)

 PlotXY(xy, clst)    

 c = BimodalCoeff(x, e)

 ChangeCell(a,b,p)

 ChangeProbe(a,b,p)

 OptionMenu(a,b, fcn, varargin)

 pos = PlaceUi(a, b ,str)

 args = RetriggerDialog(a,b, fcn)

 FitWindow(F)

 SetMenuChecks(hm, S)

 p = SelectProbe(a,b,p)

 ProbeSelector(DATA, type)

 ProbeMenu(a,b, fcn)

 res = PlotClusters(a,b,fcn)

 SetCellFromLine(a,b, cluster, cell)

 DATA = PlotCellList(DATA, varargin)

 DATA = SetCellEntry(DATA, C,  e, p, c, cellid, varargin)    

 SaveCellList(DATA)

 AddSelectorContextMenu(DATA, ax, probe)

 DATA = AddAxisContextMenu(DATA, ax)

 cmenu = AddLineContextMenu(DATA, h)

 M = CalcDistanceMatrices(DATA, nc, varargin)

 cid = AssignCluster(DATA, G)

 SetPlot(a,b, fcn)

 GMMmenu(a,b, fcn)

 GMMButton(DATA, E, fcn)

 DATA = RunGMMFit(DATA,C, varargin)

 varargout = PCCluster(a,b, fcn, varargin)

 C = GetSubCluster(Cluster, c)

 DATA = CheckTemplates(DATA, C)

 [Expt, matfile] = LoadExpt(DATA, ei)

 Expt = LoadExptA(DATA, exfile, ei)

 res = FitGaussMeans(X,N, varargin)

 [d, details]  = gmdprime(G, varargin)

 [x,y] = GetClusterXYData(DATA, p, varargin)

 x = Rprime(r)

 Cut = PlotHistogram(DATA, E, varargin)

 ApplyLayout(DATA,varargin)

 [F, isnew] = SetFigure(lb, varargin)

 XYplot(a,b,tag,fcn)

 [x, xid] = GetValues(DATA, name, dim)

 xy = PlotOneXY(DATA,names, varargin) 

 AddParameterMenu(F, callback, tag)

 CompareMean(a,b, p)

 DATA = LoadCellFile(DATA)

 [true, cellid] = isacell(DATA, ei, p)

 AddCellMenu(DATA)

 exitallv(src, evnt)

 vpts = SetVsamples(DATA, probe, np, nv)

 SetADC(a,b,fcn)

 ShowADCPos(src, data, type)

 CalcXcorr(a,b,fcn)

 FullVKeyPressed(src, ks)

 PCKeyPressed(src, ks)

 SelectTrial(src, b)

 c = TrialMarkChar(T)

 PlotOneTrial(DATA,id)

 idlist = SetTrialList(DATA)

 HistKeyPressed(src, ks)

 pos = RotateLine(pos, da)

 KeyPressed(src, ks)

 PlotGridSpikes(DATA, nspk, varargin)

 PlotQuickSpikes(DATA, nspk, varargin)

 [nt, spklst] = PlotTrialSpikes(DATA, nt, varargin)

 PlotCluster(a,b, mode)

 HistMenu(a,b, mode)

 ExptFigMenu(a,b, mode)

 PlotXcorr(a,b, pa, pb)

 Q = Cell2Cluster(cell,Clusters)

 SyncSpikes(DATA, cells)

 SpikeDraw(a,b, mode)

 xc = ShapeCorr(P,Q)

 cells = PlotAllXCorr(DATA, DataClusters, cells, varargin)

 ReplotXcorrs(a,b, type)

 PlotClusterXY(DATA, C)

 h = AddCellLabel(DATA,e,p)

 PlotAllProbes(DATA,type, varargin)

 c = MarkAxes(ax, mark)

 SummaryHit(a,b, p)

 PlotAllMeans(DATA)

 PlotMeanSpikes(C, p, cluster, varargin)

 bad = BadCluster(C)

 DATA = UseAllEvents(DATA)

 HitXYPlot(a,b, p)

 SpikeButtonPressed(a,b)

 HitImage(a,b,p)

 SpoolAllSpikes(DATA, varargin)

 PlotProbeSpikes(DATA, Vall, p, spklist,npts,offset)

 SpoolSpikes(DATA, varargin)

 csd = GetCSD(DATA, ndiff)

      DATA = SetDefaults(DATA, defmode, varargin) 
             ScrollV(src, evnt)
             ScrollSpikes(src, evnt)

 PlotFeatures(DATA, a, b, type, id, colors, C, varargin)        

 PlotPCs(pcs, a,b, type, id, colors, C,varargin)

 r = CalcRadius(E,xy)

 PlotVals(DATA, a,b, type, id, colors, varargin)

 DATA = ClassifyAll(DATA, force,varargin)

 need = NeedClusterData(Cluster, ci)

 D = CondenseCluster(C)

 DATA = IterateFit(DATA, niter)

 DATA = SetTemplateData(DATA,  cl, varargin)

 DATA = ClassifyAndFit(DATA)

 cluster = ClassifyFit(DATA, E, cnum)

 r = CalcClusterDistance(cluster, xy)

 X = GetDataStruct(DATA, f)

 [cl, cluster, xy] = ClassifySpikes(DATA, E, varargin)

 QuickSpks(DATA,nspk)

 cluster = PlotTriggerHist(DATA, cluster,varargin)

 [chspk, csdspk] = UseProbeList(DATA, nprobes)

 [MeanSpike, details] = PlotMeanSpike(DATA, varargin)

 [Scores, T, details] = CalcScores(DATA, MeanSpike)

 Labels = PCLabels(DATA, usestd)

 Labels = TemplateLabels(DATA, usestd)

 [out, TemplateUsed, DprimeUsed] = TemplatePlot(DATA, varargin)

 [bs, as] = CalculateTemplateDips(DATA)

 DATA = ReplotPCs(DATA,E, varargin)

 clusterplot = GetClusterPlots(DATA,E, plots,pt)

 h  = SetClusterIcon(DATA)

 clusterplot = GetClusterPlot(DATA,E,plots, pt)

 AddMarkToPCs(pos, space, plots, varargin)

 PlotVarE(DATA)

 DATA = RestrictTimeRange(DATA, t)

 [DATA, AllV] = ExcludeTrials(DATA, trials, add)

 SetCellCompare(a,b, cellid)

 AddCellMean(DATA, cellid)

 FinishSpikePlot(DATA)

 myAllV = GetAllV(DATA)

 ph = PlotSpikes(DATA,spkid, varargin)

 ShowFullV(src,b, fcn)

 [Vall, id] = PlotFullV(DATA, t, varargin)

 OldSetMenuCheck(F, tag, value)

 SetGUI(DATA, varargin)

 PlotTemplateScores(DATA, TemplateScores, probes)

 [dp, res] = MaxDprime(x, varargin)

  sgn = CheckSign(C, x, energy)

 [dip, details] = oldFindDip(values, energy, varargin)

 dp = CalcDprime(x, y)

 cname = ClusterFile(name,varargin)

 HistButtonPressed(src, data)

 HistButtonReleased(src, data)

 C = ClusterInfo(DATA)

 HistButtonDragged(src, data)

 h= oldDrawEllipse(E,varargin)

 h= oldDrawLine(E,varargin)

 [h, details] = DrawEllipse(E,varargin)

 h= DrawLine(E,varargin)

 [distance, cluster] = FindNearestCluster(DATA, pos)

 ButtonPressed(src, data)

 len = LineLength(l)

 distance = DistanceToEllipse(E, pos, varargin);

 in = InGraph(pt, ax)

 ButtonReleased(src, data)

 ScrollWheel(src, evnt)

 DATA = RotateCluster(DATA, angle)

 PCButtonDragged(src, data)

 pos = xyr2pos(xyr, aspect)

 C = GetClusterDef(cluster, n, varargin)

 DATA = SetEllipseDrawing(DATA, shape,varargin)

 MiscMenu(a, b, type, varargin)

 SetOption(a, b, type)

 PreferencesPopup(a,b, fcn)
 
 PopupWindow(DATA, tag, callback, varargin);

 value = Text2Val(F, tag)

 PlotResult(a, b, type)

 HitTrial(data,b, cell)

 [Expt, plotres] = PlotExptCounts(DATA)

 PlotMenu(a, b, type, varargin)

 sdx = PlotXYSequence(DATA, probe, varargin)

 PlotISI(a, b, type)

 HitISI(a,b, t)

 [isis, trials, spkids] = CalcISI(Trials, varargin)

 [C, fits] = OptimizeEllipse(DATA)

 [SSD,dipr, details ] = EllipseDip(params, DATA, state)

 guess = MinimiseEllipseDip(C, DATA)

 [SSD, details ] = MinimiseEllipse(params, DATA, state)

 [SSD, details ] = MinimiseEllipseb(params, DATA, state)

 [C, fits] = OptimizeLine(DATA)

 Plot2GaussFit(params, DATA, state)

 [SSD, details ] = MinimiseLine(params, DATA, state)
DATA = SetSpaceFromAxis(DATA, ax, varargin)

 PlotGauss2Fit(fits);

 [dp, fits, details] = Fit2Gauss(C, r, DATA, varargin)

 DATA = LoadComments(DATA)

 AddComment(a,b,str)

 C= CondenseClusters(C, go, varargin)
[AllVoltages, details] = Spikes2V(S, Cluster, allnprobes, varargin)
 p = GetProbeFromName(name)
   V = Spk2FullV(DATA, t, varargin)
       ShowTaggedProbes(DATA)

 DATA = QuickAutoCut(a,b, varargin)
    M = CalcMyDistanceMatrix(DATA, Cluster, varargin)

    end
end
