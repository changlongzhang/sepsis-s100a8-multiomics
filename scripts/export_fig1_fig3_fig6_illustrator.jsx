#target illustrator
(function(){
  app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
  var b="D:/桌面/sepsis/文章/返修投稿新主图/";
  function run(n){
    var pdf=new File(b+n+".pdf");var d=app.open(pdf);
    var t=new ExportOptionsTIFF();t.imageColorSpace=ImageColorSpace.RGB;
    t.resolution=300;t.antiAliasing=AntiAliasingMethod.ARTOPTIMIZED;
    t.lZWCompression=true;t.byteOrder=TIFFByteOrder.IBMPC;t.embedICCProfile=true;
    t.saveMultipleArtboards=true;t.artboardRange="1";
    var tf=new File(b+n+".tif");d.exportFile(tf,ExportType.TIFF,t);
    var alt=new File(b+n+"-01.tif");if(alt.exists){if(tf.exists)tf.remove();alt.rename(n+".tif");}
    var ao=new IllustratorSaveOptions();ao.pdfCompatible=true;ao.compressed=true;
    d.saveAs(new File(b+n+".ai"),ao);d.close(SaveOptions.DONOTSAVECHANGES);
    var tmp=[new File(b+n+".png"),new File(b+n+".svg")];
    for(var i=0;i<tmp.length;i++)if(tmp[i].exists)tmp[i].remove();
  }
  run("Fig1");run("Fig3");run("Fig6");
  var junk=[new File(b+"BD9317A.png"),new File(b+"BD9317E.png")];
  for(var j=0;j<junk.length;j++)if(junk[j].exists)junk[j].remove();
  app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
