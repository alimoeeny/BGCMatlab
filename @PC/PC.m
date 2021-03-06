classdef PC
properties (Constant)
CurrentVersion = '1'
STIM_GRATING = 3;
MAXCELLS = 99;
end
methods (Static)        
                     cmenu = AddCellContextMenu(DATA, type)
                   handles = AddCellLabels(DATA, eid, probe, varargin)%
                           ch = AddClusterString(DATA, h, cl)%
                                  AddComment(a,b,str)%
                     cmenu = AddContextMenu(DATA, type,varargin)%
                         C = AddCloseMenu(F, DATA, varargin);
                      DATA = AddError(DATA, varargin)%
                                  AddExptList(hm, callbacklabel, DATA)%
                             C = AddFits(DATA, C);
                       DATA = AddKey(DATA, tag, varargin)%
                       cmenu = AddLineContextMenu(DATA, h, e, p, varargin)%
                                  AddPlotMenu(sm, type)%
                       DATA = AddSelectedCells(DATA)%
                  [S, Spks] = AddSpikes(DATA, S, e,p, varargin)%
                             h = AddTitle(DATA, C, titlemode, varargin)
                                  AddTrigHist(DATA, C, cl, varargin)
                       DATA = ApplyConfig(DATA, varargin)
                       DATA = ApplyLayout(DATA, varargin)
                       DATA = AutoCompare(DATA, varargin)
                           xy = AxPos(ax, pos)
     [theta, c, details] = BestAngleGM(xy, G, dipfit, varargin)
                           im = BuildShapeMap(DATA, C, Clusters, celln, varargin)
                       DATA = CalcCellMeans(DATA)
                       DATA = CalcDistances(DATA)
                         [r,y] = CalcRadius(E,xy)
                           xc = CalcTemplateXcorr(tmpl, C)
                          res = CallAllVPcs(DATA, eid, pid, varargin)
                             Q = Cell2Cluster(cell,Clusters)
                                  CellChanged(DATA)
                                  CellFinder(DATA)
                                  CellFinderTest(DATA, type)
                                  CellTrackMenu(a,b, type)
                                  ChangeTag(a,b,fcn)
                        DATA = Check(DATA, type, varargin);
                             R = CheckAllRateSequences(DATA)
                   [X, checks] = CheckAlignment(DATA, finddata, varargin)
                    [S, new] = CheckAllSpikes(DATA, e,p, varargin)
                             F =  CheckAutoList(DATA, varargin)
                       DATA = CheckCellList(DATA,eid,p, varargin)
                       DATA = CheckClusterLineSign(DATA)
    [Clusters, DATA, e] = CheckClusterLoaded(DATA,e, varargin)
              [X, varargout] = CheckClusters(DATA, type)
                       DATA = CheckExclusion(DATA, varargin)
             [err, counts] = CheckExptRateSequence(Expt)
                       DATA = CheckExpts(DATA, type)
                          res = CheckExptsRateSequence(Expts)
                          bad = CheckFitDim(C)
                                  CheckRates(DATA, expname, cell)
                           go = CheckSpikeFile(DATA, C)
                          res = CheckSpikeFiles(DATA, type)
                             C = ClassifySpikes(DATA, E, mode)
                       DATA = ClearSelections(DATA, force, setcurrent)
                                  CloseTag(tag)
                       cells = ClusterList(DATA, Clusters)
                        need = ClusterNeedsRefresh(C, I)
                        Expt = combine(DATA, varargin);
                       DATA = CompareAutoClusters(DATA, expts)
                       DATA = CompareQuickSpks(DATA, varargin)
                                  CompareCells(DATA, Clusters,C, cella, cellb)
                 varargout = CompareCellsShape(DATA, CellList, Clusters, celid, varargin) 
                             X = CompareClusters(Ca, Cb, DATA, varargin)
                                  CompareProbes(DATA, type)
                     result = CompareProbesShape(DATA, ex, type)  
                                  CompareShape(DATA, Clusters, spk)
                       DATA = CompareShapes(DATA, type)    
                 varargout = CompareSpikeShape(A,varargin)
                    [CC, GM] = CondenseCluster(C, varargin)
                    [CC, GM] = CondenseClusters(C, varargin)
                       DATA = ConditionalPlotXY(DATA, C, force, varargin)     
                       DATA = ConvertExclusion(DATA)
                        Expt = CountExptSpikes(DATA,Expt,C,clnum, varargin)
                             E = CountSpikes(Expt, C, clnum,xcl,varargin)
                         [e,p] = cPoint(DATA)
                             e = CurrentExpt(DATA);
                     [nx, ny] = Data2Norm(x,y)
                       DATA = DeleteCell(DATA, e, p, cl)
                                  DeleteCellFromLine(a,b, cluster, cell)
                [d, drop, C] = DistanceMeasure(C, cl, type, varargin)
                  distance = DistanceToEllipse(E, pos);
                          out = doentry(DATA, Clusters, id)
                             h = DrawBox(ex, p, imtype, varargin)
                             h = DrawBoxes(DATA, imtype, varargin)
                              h= DrawEllipse(E,varargin)
                              h= DrawLine(E,varargin)
                                  EstimateDrift(DATA, where)
                       DATA = ExcludeTrialsForCell(DATA, probe, cluster, varargin)
                             s = ExLabel(DATA, j, varargin)
                                  FastAxes(ax)
                        name = FileName(DATA, ex, probe, type)
                       DATA =  FillCellList(DATA, mode)
          [eid, cid, clid] = FindCell(DATA, cellid, expts)
         [DATA, details] = FindCells(DATA, varargin);
                                  FindDuplicates(DATA, cell)
               [x, details] = FindExcludedTrials(DATA,e,p, cluster, C)
                       DATA = FindMissing(DATA)
                       shakes = FindShake(FX, varargin)
                           id = FindSpikes(DATA,C, xcl)
                  [X, xid, D] = FindXCorr(xcorrs, id)                           
                                  FinishXYPlot(ax, DATA, e,p)
     [dp, fits, details] = Fit2Gauss(C, r, varargin)
                                  FitDriftMeanSpike(M, varargin)
                                  FixBoundary(DATA, Clusters, e, p)
                            dup = FixDoubleCells(DATA, varargin)
                  Clusters = FixClusters(Clusters)
                       DATA = FlipLineCrit(a,b)
                                  Get2DMaxima(DATA)
                [value, it] = GetCheck(tag, varargin)
                      [C, ok] = GetClusterInfo(C, cl, varargin)
                                  GetComment(DATA,a,b)
               [cellid, it] = GetCurrentCell(DATA,varargin)
                             C = GetCurrentCluster(DATA)
                          eid = GetExptid(DATA, e, type, varargin);
                             e = GetExptno(DATA, eid, type, varargin);
                          tag = GetFigureTag(src)
                       mspid = GetMeanSpikeProbe(C, p)
          [str, value, it] = GetPopString(tag, varargin)
                    [value] = GetValue(DATA,type,varargin);
             [yl, details] = GetYRange(h, xl)
                                  GuiMenu(a,b, type, varargin)
                                  HitExptPlot(src, b, type, e)
                                  HitExptPoint(a,b, ex, cell)
                                  HitImage(src,b, type, varargin)
                                  HitPoint(a,b, C, mode, pt)
                                  HitPopPoint(a,b, ex, p, cell, clid)
                                  HitShapePlot(a,b,id, c)
                                  HitTrial(data,b, cell)
                                  HitXcorr(a,b, id, ex, cells)
                                  HitXcorrAll(a,b, type, cells, expts)
                                  HitXYPlot(src, b, e,p)
                                  HitXYseq(a,b)
                           in = InGraph(pt, ax)
                       DATA = InitInterface(DATA)
    [true, cellid, clid] = isacell(DATA, row, p, clid, varargin)
                          dup = isduplicate(DATA, row, p, cl, varargin)
                [same, score] = IsSameCell(xc, varargin);
                                  KeyPressed(src, ks, fcn)
                                  KeyReleased(src, ks, fcn)
                       DATA = LoadAll(a,b, type, varargin)
                                  LoadAllSpikes(DATA, varargin)
                       DATA = LoadAuto(DATA, varargin);
                       DATA = LoadCellFile(DATA, varargin)
    [Clusters, details] = LoadClusters(name)
            [Expt, DATA]  = LoadExpt(DATA, e)    
           [DATA, Expts]  = LoadExpts(DATA, varargin)
                       DATA = LoadExtra(DATA, force)
                       DATA = LoadFullVs(a,b, type, varargin)
                                  LoadSelectedSpikes(DATA,eid,pid)
                       DATA = LoadTrialLists(DATA)
                       DATA = LoadXcorrFiles(DATA)    
       [CellId, details] = MakeCellId(cellps)
       
                       pos = MakeTriplet(X, varargin)
                 varargout = MakeVarMatrix(DATA, Clusters, varargin)
                       DATA = MarkCurrentCluster(DATA)
                       h = MarkSelections(DATA, selected, varargin)
                       MarkExpts(DATA,type)
                                  MarkTrialStarts(Expt, ticks, xcl)
       [xc, synci, allxc] = meanccf(DATA, id, a,b)
                             d = meanmahal(DATA, id, a)
                           ms = MeanSpike(DATA, id, a)
                     missed = MissedCell(DATA, pos)    
                     [ex, p, cl] = Mouse2Expt(src, data);
                                  NewCellSummary(a,b,e,p)
                       DATA = NextButton(mn,b, dir)
                             s = OldExLabel(DATA, j)
                 varargout = OptionMenu(a, b, fcn, varargin)
                 varargout = PlotAllCell(DATA, type, varargin)
                                  PlotAllCellMean(DATA, type, varargin)
                                  PlotAllCellSpikes(DATA)
                                  PlotAllCellXY(DATA, varargin)
                                  PlotAllClusterMean(DATA, type, varargin)
                       DATA = PlotAllClusters(mn,b, varargin)
                                  PlotAllExpt(DATA, type)
                                  PlotAllExptProbeMean(DATA, type, varargin)
                       DATA = PlotAllProbe(DATA, type)
                                  PlotAllProbeMean(DATA, type, varargin)
                       DATA = PlotAllProbeXY(DATA,varargin)
                        args = PlotArgs(DATA)
                       DATA = PlotCellList(DATA, varargin)
                                  PlotCellRates(DATA,type)
                                  PlotCellSummary(DATA, e, p)
                     result = PlotClusterHistogram(DATA, C, refit, varargin)
                       plots = PlotClusterPoints(C, uid, cid, varargin)
         [check, Expts]  = PlotClusterRates(DATA, type,varargin)
                 varargout = PlotClusters(name, varargin)  
                       plots = PlotClusterXY(DATA, C, varargin)
                        Expt = PlotCombinedExpt(DATA, varargin)
                                  PlotCorrelogram(C, varargin)
                                  PlotDuplicates(DATA, varargin);
                                  PlotExptCells(DATA, type)
                        Expt = PlotExptCounts(DATA, e, p, cl, varargin)
                                  PlotExpts(DATA)
                       DATA = PlotExptsProbe(DATA, type);
             [x, nsp, crit] = PlotHist(xy, varargin)
                                  PlotISIHist(DATA, c);
                                  PlotMahalImage(DATA, type, varargin)
                             h = PlotMeanSpike(C, p, cluster, varargin)
                                  PlotMenu(a, b, fcn, type, varargin)
      [handles, details] = PlotPopPoints(X,Y, varargin)
                             h = PlotProbeMeans(C,type, varargin)
                                  PlotRateCheck(R)
                          Expts = PlotRateSequence(DATA, pos, varargin)
         [Expt, AllExpt] = PlotSelectedExpts(DATA, varargin)
                                  PlotShape(DATA, spk)
                       DATA = PlotShapes(DATA, type)
                       ploth = PlotSpikes(DATA, pos, spkid, Spks, C, varargin)
                                  PlotSyncSpikes(DATA, eid, probes, clnum, varargin)
                                  PlotTimeRanges(DATA)
                           nt = PlotTrialSpikes(DATA, nt, varargin)
                                  PlotXcorr(a,b, pa, pb)
                                  PlotXcorrs(DATA, xcorrs, expts, bycell)
                          sdx = PlotXYSequence(DATA, pos, varargin)
                             X = PopupWindow(a, b, fcn, varargin)
                                  PrintComments(DATA,e,p)
                          str = ProbeLabel(p, DATA, cl)
           [h, details, Spks] = QuickSpikes(DATA, pos, varargin)
                                  RateMenu(src, ks, fcn)
                                  RatePlotMenu(a,b,fcn)
                                  RateSeqKeyPressed(src, ks, fcn)
                       DATA = RateZoom(DATA,inout, cell)
                       rawxy = RawXY(C, xy)
       [DATA, Clusters] = ReadClusterResults(DATA, Clusters)
                       DATA = ReadClustersIntoData(DATA, Clusters, exid)
                       DATA = ReadTemplateResults(DATA, nc)
                                  ReFit1D(DATA, ratio)
                       DATA = ReFit3means(DATA, ratio)
                       DATA = ReFitAll(DATA, fittype)
                       DATA = ReFitGMDip(DATA, varargin)
                  Clusters = ReloadClusters(DATA, eid)
                          out = ReplotXcorrs(a,b, type, varargin)
                              ReplotXY(a,b,eid, probe,cid)
                          x = Rprime(r)
                       DATA = SaveCellList(DATA, varargin)
                       DATA = SaveCluster(DATA, pt, quick, varargin)
                                  SaveExtras(DATA)
                                  SavePlotClustersConfig(DATA,file, varargin);
                                  ScrollWheel(src, evnt)
         [cells,e,p,clid] = SelectedCells(selectid, CellList, varargin);
                                  SelectTrial(src, b, type)
                       DATA = SetAllXYClustering(DATA, onoff) 
                       DATA = SetCellEntry(DATA, C,  e, p, c, cellid, varargin)    
                                  SetCellFromLine(a,b, cluster, cell, varargin)
                                  SetCellFromSubplot(a,b, cell)
                                  SetCellNumber(a,b, fcn)
                                  SetCheckExclusive(a)
                       DATA = SetDefaults(DATA)
                       DATA = SetExpt(mn,b,  fcn)
            [DATA, Expts] = SetExptList(DATA, Expts, varargin)
                       Expts = SetExptTimeOffset(Expts)
                  [f, isnew] = SetFigure(DATA, tag, varargin)
                                  SetGUI(DATA)
                    plotpos = SetPlotPos(DATA,np, nr, nc)
                       DATA = SetProbe(mn,b, dir)
                                  SetRateseqPlot(a,b,cellid)
                                  SetTrialList(DATA, C, strial)
          varargout = SetValue(DATA, type, value, varargin);
          [voffset, ylim] = SetVOffset(DATA, AllSpikes, e, p)
             [xc, details] = ShapeCorr(P,Q, varargin)
                          xcs = ShiftXcorr(allshape, a, b, npts)
                      DATA  = ShowData(DATA, ex,p, varargin)
     [sizes, d, quality] = SpikeDistance(DATA, eid)
                                  SpikeDistances(DATA, type)
                        name = SpikeFileName(DATA, C)
                    stopped = SpoolAllProbes(DATA, e, Spikes, Clusters, varargin)
                       DATA = SpoolCurrentSpikes(mn,b, varargin)
                          res = SpoolSpikeFile(DATA,e,p)
               [stopped, h] = SpoolSpikes(DATA, pos, varargin)
                       synci = SyncIndices(xc)
                                  TagMenu(a, b, fcn)
                              h= testDrawEllipse(E,varargin)
                              h= testDrawLine(E,varargin)
                                  TightPlot(ax)
                       DATA = TrackCluster(DATA, spk, ex, varargin)
                                  Update(a,b)
                             C = UpdateWithNewCut(C, DATA)
                             WindowMotion(src, data, type)
                                  XYButtonDragged(src, data)
                                  XYButtonPressed(src, data)
                                  XYButtonReleased(src, data)
                                  XYCluster(src,b, type, varargin)
                                  XYKeyPressed(src, ks, fcn)
                           xy = XYSpace(C)
end
end
