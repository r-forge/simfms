<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

<head>
  <title>simfms: Simulation of clustered multi-state data</title>

  <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
  <meta name="author" content="Federico Rotolo" />
  <meta name="description" content="Simulation of clustered multi-state data" />
  <meta name="keywords" content="simulation, frailty, multi-state, copula, survival, R, package" /> 

  <link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  <link rel="stylesheet" type="text/css" href="https://r-forge.r-project.org/themes/css/gforge.css" />
  <link rel="stylesheet" type="text/css" href="https://r-forge.r-project.org/themes/rforge/css/theme.css" />
</head>
<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- project title  -->

<h2>Simulation of clustered multi-state data</h2>

<div style="width:425px" id="__ss_9422763">
  <strong style="display:block;margin:12px 0 4px">
    <a href="http://www.slideshare.net/federicorotolo/a-copulabased-simulation-method-for-clustered-multistate-survival-data" 
       title="A copula-based Simulation Method for Clustered Multi-State Survival Data"
       target="_blank">A copula-based Simulation Method for Clustered Multi-State Survival Data</a>
   </strong>
   <iframe src="http://www.slideshare.net/slideshow/embed_code/9422763" width="425" height="355" frameborder="0" marginwidth="0" marginheight="0" scrolling="no"></iframe>
   <!--div style="padding:5px 0 12px"> View more <a href="http://www.slideshare.net/" target="_blank">presentations</a>
    from <a href="http://www.slideshare.net/federicorotolo" target="_blank">federicorotolo</a>
   </div-->
</div>


<p><strong>Project summary</strong>:
  <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/">here</a>. </p>


</body>
</html>
