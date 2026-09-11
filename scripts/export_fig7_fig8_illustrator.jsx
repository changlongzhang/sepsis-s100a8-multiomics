#target illustrator
(function(){
 function run(n){
  var base="D:/桌面/sepsis/文章/返修投稿新主图/"; var master=new File(base+n+"_python_master.pdf"); var doc=app.open(master);
  var arial=null;try{arial=app.textFonts.getByName("ArialMT");}catch(e1){try{arial=app.textFonts.getByName("Arial");}catch(e2){}}
  if(arial!==null){for(var i=0;i<doc.textFrames.length;i++)doc.textFrames[i].textRange.characterAttributes.textFont=arial;}
  var t=new ExportOptionsTIFF();t.imageColorSpace=ImageColorSpace.RGB;t.resolution=300;t.antiAliasing=AntiAliasingMethod.ARTOPTIMIZED;t.lZWCompression=true;t.byteOrder=TIFFByteOrder.IBMPC;t.embedICCProfile=true;t.saveMultipleArtboards=true;t.artboardRange="1";
  var tf=new File(base+n+".tif");doc.exportFile(tf,ExportType.TIFF,t);var ab=new File(base+n+"-01.tif");if(ab.exists){if(tf.exists)tf.remove();ab.rename(n+".tif");}
  var po=new PDFSaveOptions();po.compatibility=PDFCompatibility.ACROBAT7;po.preserveEditability=true;po.optimization=true;po.viewAfterSaving=false;doc.saveAs(new File(base+n+".pdf"),po);
  var ao=new IllustratorSaveOptions();ao.pdfCompatible=true;ao.embedICCProfile=true;ao.compressed=true;doc.saveAs(new File(base+n+".ai"),ao);doc.close(SaveOptions.DONOTSAVECHANGES);
  var tmp=[master,new File(base+n+".png"),new File(base+n+".svg"),new File(base+n+"_preview.png")];for(var j=0;j<tmp.length;j++)if(tmp[j].exists)tmp[j].remove();
 }
 app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;run("Fig7");run("Fig8");app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
