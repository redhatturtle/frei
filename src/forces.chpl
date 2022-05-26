module Forces
{
  proc cp(famlIdx : int) : [] real
  {
    use Input;

    // Loop though all boundaries
    forall thisBoco in frMesh.bocolist
    {
      ref famlIdx = thisBoco.family;
      ref famlName = frMesh.famlList[famlidx].name;

      // Check this boundary belongs to one of the families we want to calculate Cp over
      for forcesFamlIdx in Input.focesFaml.domain
      {
        if famlName == Input.forcesFaml[forcesFamlIdx]
        {
          // It does
          ref faceIdx = thisBoco.face;

          for fpIdx in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
          {
            var avgsol = (frMesh.solFP[fpIdx, 1, ..] + frMesh.solFP[fpIdx, 1, ..])/2.0;
            var cp = (pres_cv(avgSol) - Input.presInf)/(Input.densInf*Input.velInf/2.0);
          }
        }
      }
    }

  return ???;
  }

  proc lift() : real
  {}

  proc drag() : real
  {}
}
