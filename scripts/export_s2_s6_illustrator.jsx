#target illustrator
(function(){
  app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
  var b="D:/桌面/sepsis/文章/返修投稿新主图/";
  function run(n){
    var m=new File(b+n+"_fig_python_master.pdf");var d=app.open(m);
    var t=new ExportOptionsTIFF();t.imageColorSpace=ImageColorSpace.RGB;
    t.resolution=300;t.antiAliasing=AntiAliasingMethod.ARTOPTIMIZED;
    t.lZWCompression=true;t.byteOrder=TIFFByteOrder.IBMPC;t.embedICCProfile=true;
    t.saveMultipleArtboards=true;t.artboardRange="1";
    var tf=new File(b+n+"_fig.tif");d.exportFile(tf,ExportType.TIFF,t);
    var alt=new File(b+n+"_fig-01.tif");
    if(alt.exists){if(tf.exists)tf.remove();alt.rename(n+"_fig.tif");}
    var po=new PDFSaveOptions();po.compatibility=PDFCompatibility.ACROBAT7;
    po.preserveEditability=true;po.optimization=true;
    d.saveAs(new File(b+n+"_fig.pdf"),po);
    var ao=new IllustratorSaveOptions();ao.pdfCompatible=true;ao.compressed=true;
    d.saveAs(new File(b+n+"_fig.ai"),ao);
    d.close(SaveOptions.DONOTSAVECHANGES);
    var tmp=[m,new File(b+n+"_fig.png"),new File(b+n+"_fig.svg")];
    for(var i=0;i<tmp.length;i++)if(tmp[i].exists)tmp[i].remove();
  }
  for(var k=2;k<=6;k++)run("S"+k);
  app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
