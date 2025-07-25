
classdef DRS_DAE3 < DAC_ODEclass
    properties (SetAccess = protected, Hidden = true)
        g  ;
        m1 ;
        m3 ;
        psi;
        %LF;
        anc;
        LFb;
        LP ;
        %Iz1;
        Iz3;
    end

    methods   
        
        function self = DRS_DAE3(g,m1,m3,psi,anc,LFb,LP,Iz3)
            neq  = 7;   %number of equations
            ninv = 0;  %number of hidden constraints

            self@DAC_ODEclass('DRS_DAE',neq,ninv);

            %setup of the parameter of the ODE
            self.g= g;
            self.m1 = m1;
            self.m3 = m3;
            self.psi = psi;
            %LF  = self.LF;
            self.anc = anc;
            self.LFb = LFb;
            self.LP = LP;
            %Iz1 = self.Iz1;
            self.Iz3 = Iz3;
        end
        function res__f = f( self, T, vars__ )

    g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

      % evaluate function
      res__1 = s__dot;
      t1 = cos(theta3);
      t5 = cos(theta1);
      res__2 = 0.1e1 / t5 / LP * theta3__dot * anc * t1;
      res__3 = 0.1e1 * theta3__dot;
      t10 = T ^ 2;
      t12 = t10 * T;
      t14 = t10 ^ 2;
      t16 = t14 * T;
      t18 = t14 * t10;
      t20 = t14 * t12;
      t22 = t14 ^ 2;
      t28 = -0.751955623237439977e4 + T * (-0.842440941625551932e4) + t10 * 0.383430218313185498e7 + t12 * (-0.792218461538717270e9) + t14 * 0.875421336981854248e11 + t16 * (-0.668171235384638281e13) + t18 * 0.362875310184053500e15 + t20 * (-0.14533139837090800e17) + t22 * 0.445733173550047232e18 + t22 * T * (-0.106507163715741614e20) + t22 * t10 * 0.199409962072002265e21;
      t39 = t22 ^ 2;
      t50 = t22 * t12 * (-0.296788012772351384e22) + t22 * t14 * 0.350712552560261826e23 + t22 * t16 * (-0.331015879728087465e24) + t22 * t18 * 0.246747972057886706e25 + t22 * t20 * (-0.143677886346597917e26) + t39 * 0.638797917538537755e26 + t39 * T * (-0.210948397437314510e27) + t39 * t10 * 0.482444576644185378e27 + t39 * t12 * (-0.688070947738516329e27) + t39 * t14 * 0.459999033329155986e27 + 4 * lambda1;
      res__4 = 0.1e1 / m1 * (t28 + t50);
      t53 = theta3 ^ 2;
      t54 = t53 * theta3;
      t55 = 0.592621945e4 * t54;
      t56 = 0.5740478308e5 * t53;
      t57 = 0.182598982e6 * theta3;
      t58 = g * m3;
      t59 = 0.1e1 * t58;
      t62 = cos(psi);
      t64 = 0.474097556e4 * t54;
      t68 = LFb * (t64 + t53 * (-0.490267877100000042e5) + theta3 * 0.167958569199998077e6 - 0.190914737699991150e6);
      t69 = sin(psi);
      t81 = sin(theta3);
      t87 = LFb ^ 2;
      t88 = t87 * m3;
      res__5 = 0.1e1 / (t88 + Iz3) * (t1 * (t62 * LFb * (t55 - t56 + t57 - t59 - 0.1908391204e6) + t69 * t68 - 0.1e1 * anc * lambda2) + t81 * (t62 * t68 + t69 * LFb * (-t55 + t56 - t57 + t58 + 0.1908391204e6) + lambda1 * anc) - 0.4266877991e3 * t54 + t53 * 0.433173388900000373e4 + theta3 * (-0.145530395300000237e5) + 0.161877700100000457e5);
      t91 = 0.6096051475e5 * t53;
      t92 = 0.2071123758e6 * theta3;
      t94 = LFb * (0.2328287627e6 - t55 + t91 - t92 + t58);
      t95 = LP ^ 2;
      t96 = t69 * t95;
      t100 = t64 + t53 * (-0.445821231199999966e5) + theta3 * 0.139256177600000025e6 - 0.145264992200000037e6;
      t101 = LFb * t100;
      t102 = t62 * t95;
      t104 = lambda1 * t95;
      t105 = anc * t104;
      t106 = t102 * t101 + t96 * t94 + t105;
      t107 = t1 ^ 2;
      t110 = -LFb * t100;
      t113 = lambda2 * t95;
      t114 = anc * t113;
      t119 = theta3__dot ^ 2;
      t120 = Iz3 * t119;
      t122 = m3 * t119;
      t125 = -0.3200158492e3 * t53 + 0.2165866944e4 * theta3 - 0.25e0 * t120 - 0.3638259882e4 - 0.25e0 * t87 * t122;
      t126 = t95 * t125;
      t132 = t53 * (-0.324880041799999981e4) + theta3 * 0.109147796500000004e5 + 0.3200158492e3 * t54 - 0.121408275100000028e5;
      t133 = t95 * t132;
      t138 = 0.75e0 * t58;
      t139 = -0.1431293403e6 + 0.4444664588e4 * t54 - 0.4305358732e5 * t53 + 0.1369492365e6 * theta3 - t138;
      t140 = LFb * t139;
      t145 = 0.143186053299993364e6 - 0.355573167e4 * t54 + t53 * 0.367700907800000059e5 + theta3 * (-0.125968926899998565e6);
      t146 = LFb * t145;
      t148 = 0.75e0 * t105;
      t150 = t5 ^ 2;
      t151 = t150 ^ 2;
      t152 = t151 * t150;
      t154 = t107 * t1;
      t155 = t154 * LP;
      t156 = anc ^ 2;
      t157 = t156 * lambda1;
      t159 = 0.25e0 * t157 * t155;
      t162 = LFb * (-0.2328287627e6 + t55 - t91 + t92 - t59);
      t165 = t96 * t101 + t102 * t162 - 0.1e1 * t114;
      t166 = sin(theta1);
      t168 = t81 * LP;
      t169 = t156 * lambda2;
      t176 = -t95 * t132;
      t185 = LFb * (-0.118524389e4 * t54 + t53 * 0.781203233999999793e4 + theta3 * (-0.132872507600000008e5) + 0.207893892499999993e4);
      t192 = LFb * (0.896994224e5 + 0.1790692744e5 * t53 - 0.1481554862e4 * t54 - 0.7016313935e5 * theta3 + 0.25e0 * t58);
      t198 = t151 * t5;
      t200 = t156 * t166;
      t201 = lambda2 * t155;
      t204 = lambda1 * t168;
      t214 = -t95 * t125;
      t219 = -LFb * t139;
      t222 = -LFb * t145;
      t227 = lambda2 * t81;
      t232 = t150 * t5;
      t235 = anc * t54;
      t237 = anc * t53;
      t239 = theta3 * anc;
      t242 = LFb * (0.1431860533e6 * anc - 0.355573167e4 * t235 + 0.3677009078e5 * t237 - 0.1259689269e6 * t239);
      t243 = t69 * LP;
      t245 = 0.4444664588e4 * t235;
      t246 = 0.430535873e5 * t237;
      t247 = 0.1369492365e6 * t239;
      t248 = 0.1431293403e6 + t138;
      t251 = LFb * (anc * t248 - t245 + t246 - t247);
      t252 = t62 * LP;
      t254 = lambda2 * LP;
      t260 = 0.75e0 * t119 * t156 * Iz3;
      t263 = 0.75e0 * t122 * t156 * t87;
      t269 = LFb * (-anc * t248 + t245 - t246 + t247);
      t272 = lambda1 * LP;
      t282 = LP * (0.3200158492e3 * t235 - 0.3248800418e4 * t237 + 0.1091477965e5 * t239 - 0.1214082751e5 * anc);
      t293 = LP * (0.75e0 * anc * m3 * t87 * t119 + 0.75e0 * anc * t120);
      t295 = t1 * t166 * t81;
      t299 = -t260 - t263;
      t308 = m1 * t156;
      t312 = (0.2e1 * m1 * t107 * t156 - 0.4e1 * t88 - 0.1e1 * t308 - 0.4e1 * Iz3) * t95;
      t315 = t156 * t95;
      t316 = t1 * t315;
      t317 = t166 ^ 2;
      t324 = t166 * m1;
      res__6 = 0.4000000000e1 / (t198 * t166 * t312 + 0.2e1 * t151 * t81 * t317 * m1 * t316 - 0.1e1 * t232 * t324 * t107 * t315) * theta3__dot * anc * (t152 * (t107 * t106 + t1 * (t81 * (t102 * t94 + t96 * t110 + t114) + t126) + t81 * t133 + t96 * t140 + t102 * t146 - t148) + t198 * (-t159 + t107 * (t166 * t165 + 0.25e0 * t169 * t168) + t1 * t166 * (t81 * t106 + t176) + t166 * (t81 * t126 + t96 * t185 + t102 * t192 + 0.25e0 * t114)) + t151 * (-0.25e0 * t201 * t200 + t107 * (-0.25e0 * t204 * t200 + t96 * t162 + t102 * t110 - 0.1e1 * t105) + t1 * (t81 * t165 + t214) + t81 * t176 + t96 * t219 + t102 * t222 + t148) + t232 * (t159 - 0.25e0 * t156 * t227 * t107 * LP) + t150 * (t154 * (t166 * (t156 * t254 + t243 * t242 + t252 * t251) + t260 + t263) + t107 * t166 * (t81 * (t243 * t269 + t252 * t242 - 0.75e0 * t156 * t272) + t282) + t295 * t293) + t154 * t299) * m1;
      t332 = m1 * t95;
      t333 = t69 * t332;
      t335 = t62 * t332;
      t337 = anc * m1;
      t338 = t337 * t104;
      t339 = t335 * t101 + t333 * t94 + t338;
      t343 = t337 * t113;
      t349 = t81 * m1;
      t353 = 0.75e0 * t338;
      t362 = t333 * t101 + t335 * t162 - 0.1e1 * t343;
      t373 = t308 * t272;
      t392 = t156 * t324;
      t402 = t308 * t254;
      t427 = LP * m1;
      t428 = t69 * t427;
      t430 = t62 * t427;
      res__7 = 0.4000000000e1 / (t152 * t312 + 0.2e1 * t198 * t81 * t324 * t316 - 0.1e1 * t151 * m1 * t107 * t315) * theta3__dot * anc * (t152 * (t107 * t339 + t1 * (t81 * (t333 * t110 + t335 * t94 + t343) + m1 * t126) + t349 * t133 + t333 * t140 + t335 * t146 - t353) + t198 * (0.25e0 * m1 * t157 * t155 + t107 * (t166 * t362 - 0.25e0 * m1 * t169 * t168) + t1 * (t166 * (m1 * t176 + t81 * t339) - 0.25e0 * t373 + LP * (-0.1e1 * lambda1 * Iz3 - 0.1e1 * t87 * m3 * lambda1)) + t166 * (t349 * t126 + t333 * t185 + t335 * t192 + 0.25e0 * t343)) + t151 * (0.25e0 * t201 * t392 + t107 * (0.25e0 * t204 * t392 + t333 * t162 + t335 * t110 - 0.1e1 * t338) + t1 * (t166 * (-0.25e0 * t402 + LP * (-0.1e1 * lambda2 * Iz3 - 0.1e1 * t87 * m3 * lambda2)) + t81 * t362 + m1 * t214) + t349 * t176 + t333 * t219 + t335 * t222 + t353) + 0.25e0 * t308 * t227 * t107 * LP * t232 + t150 * (t154 * (t166 * (t428 * t242 + t430 * t251 + 0.75e0 * t402) - m1 * t299) + t107 * t166 * (t81 * (t428 * t269 + t430 * t242 - 0.75e0 * t373) + m1 * t282) + t295 * m1 * t293) + t154 * m1 * t299);
      

      % store on output
      res__f = zeros(7,1);
      res__f(1) = res__1;
      res__f(2) = res__2;
      res__f(3) = res__3;
      res__f(4) = res__4;
      res__f(5) = res__5;
      res__f(6) = res__6;
      res__f(7) = res__7;

    end
        function res__DfDx = DfDx( self, T, vars__ )

    g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;
      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

      % evaluate function
      res__1_1 = s__dot;
      t1 = cos(theta3);
      t5 = cos(theta1);
      res__2_1 = 0.1e1 / t5 / LP * theta3__dot * anc * t1;
      res__3_1 = 0.1e1 * theta3__dot;
      t10 = T ^ 2;
      t12 = t10 * T;
      t14 = t10 ^ 2;
      t16 = t14 * T;
      t18 = t14 * t10;
      t20 = t14 * t12;
      t22 = t14 ^ 2;
      t28 = -0.751955623237439977e4 + T * (-0.842440941625551932e4) + t10 * 0.383430218313185498e7 + t12 * (-0.792218461538717270e9) + t14 * 0.875421336981854248e11 + t16 * (-0.668171235384638281e13) + t18 * 0.362875310184053500e15 + t20 * (-0.14533139837090800e17) + t22 * 0.445733173550047232e18 + t22 * T * (-0.106507163715741614e20) + t22 * t10 * 0.199409962072002265e21;
      t39 = t22 ^ 2;
      t50 = t22 * t12 * (-0.296788012772351384e22) + t22 * t14 * 0.350712552560261826e23 + t22 * t16 * (-0.331015879728087465e24) + t22 * t18 * 0.246747972057886706e25 + t22 * t20 * (-0.143677886346597917e26) + t39 * 0.638797917538537755e26 + t39 * T * (-0.210948397437314510e27) + t39 * t10 * 0.482444576644185378e27 + t39 * t12 * (-0.688070947738516329e27) + t39 * t14 * 0.459999033329155986e27 + 4 * lambda1;
      res__4_1 = 0.1e1 / m1 * (t28 + t50);
      t53 = theta3 ^ 2;
      t54 = t53 * theta3;
      t55 = 0.592621945e4 * t54;
      t56 = 0.5740478308e5 * t53;
      t57 = 0.182598982e6 * theta3;
      t58 = g * m3;
      t59 = 0.1e1 * t58;
      t62 = cos(psi);
      t64 = 0.474097556e4 * t54;
      t68 = LFb * (t64 + t53 * (-0.490267877100000042e5) + theta3 * 0.167958569199998077e6 - 0.190914737699991150e6);
      t69 = sin(psi);
      t81 = sin(theta3);
      t87 = LFb ^ 2;
      t88 = t87 * m3;
      res__5_1 = 0.1e1 / (t88 + Iz3) * (t1 * (t62 * LFb * (t55 - t56 + t57 - t59 - 0.1908391204e6) + t69 * t68 - 0.1e1 * anc * lambda2) + t81 * (t62 * t68 + t69 * LFb * (-t55 + t56 - t57 + t58 + 0.1908391204e6) + lambda1 * anc) - 0.4266877991e3 * t54 + t53 * 0.433173388900000373e4 + theta3 * (-0.145530395300000237e5) + 0.161877700100000457e5);
      t91 = 0.6096051475e5 * t53;
      t92 = 0.2071123758e6 * theta3;
      t94 = LFb * (0.2328287627e6 - t55 + t91 - t92 + t58);
      t95 = LP ^ 2;
      t96 = t69 * t95;
      t100 = t64 + t53 * (-0.445821231199999966e5) + theta3 * 0.139256177600000025e6 - 0.145264992200000037e6;
      t101 = LFb * t100;
      t102 = t62 * t95;
      t104 = lambda1 * t95;
      t105 = anc * t104;
      t106 = t102 * t101 + t96 * t94 + t105;
      t107 = t1 ^ 2;
      t110 = -LFb * t100;
      t113 = lambda2 * t95;
      t114 = anc * t113;
      t119 = theta3__dot ^ 2;
      t120 = Iz3 * t119;
      t122 = m3 * t119;
      t125 = -0.3200158492e3 * t53 + 0.2165866944e4 * theta3 - 0.25e0 * t120 - 0.3638259882e4 - 0.25e0 * t87 * t122;
      t126 = t95 * t125;
      t132 = t53 * (-0.324880041799999981e4) + theta3 * 0.109147796500000004e5 + 0.3200158492e3 * t54 - 0.121408275100000028e5;
      t133 = t95 * t132;
      t138 = 0.75e0 * t58;
      t139 = -0.1431293403e6 + 0.4444664588e4 * t54 - 0.4305358732e5 * t53 + 0.1369492365e6 * theta3 - t138;
      t140 = LFb * t139;
      t145 = 0.143186053299993364e6 - 0.355573167e4 * t54 + t53 * 0.367700907800000059e5 + theta3 * (-0.125968926899998565e6);
      t146 = LFb * t145;
      t148 = 0.75e0 * t105;
      t150 = t5 ^ 2;
      t151 = t150 ^ 2;
      t152 = t151 * t150;
      t154 = t107 * t1;
      t155 = t154 * LP;
      t156 = anc ^ 2;
      t157 = t156 * lambda1;
      t159 = 0.25e0 * t157 * t155;
      t162 = LFb * (-0.2328287627e6 + t55 - t91 + t92 - t59);
      t165 = t96 * t101 + t102 * t162 - 0.1e1 * t114;
      t166 = sin(theta1);
      t168 = t81 * LP;
      t169 = t156 * lambda2;
      t176 = -t95 * t132;
      t185 = LFb * (-0.118524389e4 * t54 + t53 * 0.781203233999999793e4 + theta3 * (-0.132872507600000008e5) + 0.207893892499999993e4);
      t192 = LFb * (0.896994224e5 + 0.1790692744e5 * t53 - 0.1481554862e4 * t54 - 0.7016313935e5 * theta3 + 0.25e0 * t58);
      t198 = t151 * t5;
      t200 = t156 * t166;
      t201 = lambda2 * t155;
      t204 = lambda1 * t168;
      t214 = -t95 * t125;
      t219 = -LFb * t139;
      t222 = -LFb * t145;
      t227 = lambda2 * t81;
      t232 = t150 * t5;
      t235 = anc * t54;
      t237 = anc * t53;
      t239 = theta3 * anc;
      t242 = LFb * (0.1431860533e6 * anc - 0.355573167e4 * t235 + 0.3677009078e5 * t237 - 0.1259689269e6 * t239);
      t243 = t69 * LP;
      t245 = 0.4444664588e4 * t235;
      t246 = 0.430535873e5 * t237;
      t247 = 0.1369492365e6 * t239;
      t248 = 0.1431293403e6 + t138;
      t251 = LFb * (anc * t248 - t245 + t246 - t247);
      t252 = t62 * LP;
      t254 = lambda2 * LP;
      t260 = 0.75e0 * t119 * t156 * Iz3;
      t263 = 0.75e0 * t122 * t156 * t87;
      t269 = LFb * (-anc * t248 + t245 - t246 + t247);
      t272 = lambda1 * LP;
      t282 = LP * (0.3200158492e3 * t235 - 0.3248800418e4 * t237 + 0.1091477965e5 * t239 - 0.1214082751e5 * anc);
      t293 = LP * (0.75e0 * anc * m3 * t87 * t119 + 0.75e0 * anc * t120);
      t295 = t1 * t166 * t81;
      t299 = -t260 - t263;
      t308 = m1 * t156;
      t312 = (0.2e1 * m1 * t107 * t156 - 0.4e1 * t88 - 0.1e1 * t308 - 0.4e1 * Iz3) * t95;
      t315 = t156 * t95;
      t316 = t1 * t315;
      t317 = t166 ^ 2;
      t324 = t166 * m1;
      res__6_1 = 0.4000000000e1 / (t198 * t166 * t312 + 0.2e1 * t151 * t81 * t317 * m1 * t316 - 0.1e1 * t232 * t324 * t107 * t315) * theta3__dot * anc * (t152 * (t107 * t106 + t1 * (t81 * (t102 * t94 + t96 * t110 + t114) + t126) + t81 * t133 + t96 * t140 + t102 * t146 - t148) + t198 * (-t159 + t107 * (t166 * t165 + 0.25e0 * t169 * t168) + t1 * t166 * (t81 * t106 + t176) + t166 * (t81 * t126 + t96 * t185 + t102 * t192 + 0.25e0 * t114)) + t151 * (-0.25e0 * t201 * t200 + t107 * (-0.25e0 * t204 * t200 + t96 * t162 + t102 * t110 - 0.1e1 * t105) + t1 * (t81 * t165 + t214) + t81 * t176 + t96 * t219 + t102 * t222 + t148) + t232 * (t159 - 0.25e0 * t156 * t227 * t107 * LP) + t150 * (t154 * (t166 * (t156 * t254 + t243 * t242 + t252 * t251) + t260 + t263) + t107 * t166 * (t81 * (t243 * t269 + t252 * t242 - 0.75e0 * t156 * t272) + t282) + t295 * t293) + t154 * t299) * m1;
      t332 = m1 * t95;
      t333 = t69 * t332;
      t335 = t62 * t332;
      t337 = anc * m1;
      t338 = t337 * t104;
      t339 = t335 * t101 + t333 * t94 + t338;
      t343 = t337 * t113;
      t349 = t81 * m1;
      t353 = 0.75e0 * t338;
      t362 = t333 * t101 + t335 * t162 - 0.1e1 * t343;
      t373 = t308 * t272;
      t392 = t156 * t324;
      t402 = t308 * t254;
      t427 = LP * m1;
      t428 = t69 * t427;
      t430 = t62 * t427;
      res__7_1 = 0.4000000000e1 / (t152 * t312 + 0.2e1 * t198 * t81 * t324 * t316 - 0.1e1 * t151 * m1 * t107 * t315) * theta3__dot * anc * (t152 * (t107 * t339 + t1 * (t81 * (t333 * t110 + t335 * t94 + t343) + m1 * t126) + t349 * t133 + t333 * t140 + t335 * t146 - t353) + t198 * (0.25e0 * m1 * t157 * t155 + t107 * (t166 * t362 - 0.25e0 * m1 * t169 * t168) + t1 * (t166 * (m1 * t176 + t81 * t339) - 0.25e0 * t373 + LP * (-0.1e1 * lambda1 * Iz3 - 0.1e1 * t87 * m3 * lambda1)) + t166 * (t349 * t126 + t333 * t185 + t335 * t192 + 0.25e0 * t343)) + t151 * (0.25e0 * t201 * t392 + t107 * (0.25e0 * t204 * t392 + t333 * t162 + t335 * t110 - 0.1e1 * t338) + t1 * (t166 * (-0.25e0 * t402 + LP * (-0.1e1 * lambda2 * Iz3 - 0.1e1 * t87 * m3 * lambda2)) + t81 * t362 + m1 * t214) + t349 * t176 + t333 * t219 + t335 * t222 + t353) + 0.25e0 * t308 * t227 * t107 * LP * t232 + t150 * (t154 * (t166 * (t428 * t242 + t430 * t251 + 0.75e0 * t402) - m1 * t299) + t107 * t166 * (t81 * (t428 * t269 + t430 * t242 - 0.75e0 * t373) + m1 * t282) + t295 * m1 * t293) + t154 * m1 * t299);
      
      % store on output
      res__DfDx = zeros(7,1);
      res__DfDx(1,1) = res__1_1;
      res__DfDx(2,1) = res__2_1;
      res__DfDx(3,1) = res__3_1;
      res__DfDx(4,1) = res__4_1;
      res__DfDx(5,1) = res__5_1;
      res__DfDx(6,1) = res__6_1;
      res__DfDx(7,1) = res__7_1;
    end
        function res__h = h( self, T, vars__ )

    g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;
      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

      % evaluate function
      t1 = sin(theta1);
      t3 = sin(theta3);
      res__1 = t1 * LP - anc * t3 + y1 - y4;
      t5 = cos(theta1);
      t7 = cos(theta3);
      res__2 = t5 * LP - anc * t7 + L1 + s + x1 - x4;
      res__3 = -(lambda1 * t1 - lambda2 * t5) * LP;
      res__4 = -s__dot + theta3__dot / t5 * (t7 * t1 - t3 * t5) * anc;
      t20 = theta3 ^ 2;
      t21 = t20 * theta3;
      t22 = 0.592621945e4 * t21;
      t23 = 0.5740478308e5 * t20;
      t24 = 0.182598982e6 * theta3;
      t25 = g * m3;
      t29 = m1 * anc * (t22 - t23 + t24 - 0.1e1 * t25 - 0.1908391204e6);
      t30 = LP * LFb;
      t31 = sin(psi);
      t32 = t31 * t30;
      t37 = -0.474097556e4 * t21 + t20 * 0.490267877100000042e5 + theta3 * (-0.167958569199998077e6) + 0.190914737699991150e6;
      t39 = m1 * anc * t37;
      t40 = cos(psi);
      t41 = t40 * t30;
      t44 = anc ^ 2;
      t45 = m1 * t44;
      t48 = t32 * t29 + t41 * t39 - 0.1e1 * t45 * LP * lambda1;
      t49 = t7 ^ 2;
      t53 = -m1 * anc * t37;
      t57 = t45 * LP * lambda2;
      t61 = theta3__dot ^ 2;
      t64 = LFb ^ 2;
      t66 = t64 * m3 + Iz3;
      t74 = -0.4266877991e3 * t21 + t20 * 0.433173388900000009e4 + theta3 * (-0.145530395300000000e5) + 0.161877700100000038e5;
      t76 = LP * m1;
      t81 = m1 * anc * (-t22 + t23 - t24 + t25 + 0.1908391204e6);
      t84 = T ^ 2;
      t85 = t84 ^ 2;
      t86 = t85 ^ 2;
      t87 = t86 ^ 2;
      t88 = t87 * t85;
      t91 = t84 * T;
      t92 = t87 * t91;
      t95 = t87 * t84;
      t98 = t87 * T;
      t103 = t85 * t91;
      t104 = t86 * t103;
      t107 = t85 * t84;
      t108 = t86 * t107;
      t111 = t85 * T;
      t112 = t86 * t111;
      t115 = t86 * t85;
      t118 = t86 * t91;
      t122 = t86 * t84;
      t125 = t86 * T;
      t145 = -0.7519556232e4 + 0.4e1 * lambda1;
      t147 = 0.1994099621e21 * m3 * t122 - 0.1065071637e20 * m3 * t125 + 0.4457331736e18 * m3 * t86 - 0.1453313984e17 * m3 * t103 + 0.3628753102e15 * m3 * t107 - 0.6681712354e13 * m3 * t111 + 0.8754213370e11 * m3 * t85 - 0.7922184615e9 * m3 * t91 + 0.3834302183e7 * m3 * t84 - 0.8424409416e4 * T * m3 + m3 * t145;
      t170 = t64 * (0.4599990333e27 * m3 * t88 - 0.6880709477e27 * m3 * t92 + 0.4824445766e27 * m3 * t95 - 0.2109483974e27 * m3 * t98 + 0.6387979175e26 * t87 * m3 - 0.1436778863e26 * m3 * t104 + 0.2467479721e25 * m3 * t108 - 0.3310158797e24 * m3 * t112 + 0.3507125526e23 * m3 * t115 - 0.2967880128e22 * m3 * t118 + t147) + m1 * t44 * lambda1 + 0.4599990333e27 * Iz3 * t88 - 0.6880709477e27 * Iz3 * t92 + 0.4824445766e27 * Iz3 * t95 - 0.2109483974e27 * Iz3 * t98 + 0.6387979175e26 * Iz3 * t87 - 0.1436778863e26 * Iz3 * t104 + 0.2467479721e25 * Iz3 * t108 - 0.3310158797e24 * Iz3 * t112 + 0.3507125526e23 * Iz3 * t115;
      t194 = -0.2967880128e22 * Iz3 * t118 + 0.1994099621e21 * Iz3 * t122 - 0.1065071637e20 * Iz3 * t125 + 0.4457331736e18 * Iz3 * t86 - 0.1453313984e17 * Iz3 * t103 + 0.3628753102e15 * Iz3 * t107 - 0.6681712354e13 * Iz3 * t111 + 0.8754213370e11 * Iz3 * t85 - 0.7922184615e9 * Iz3 * t91 + 0.3834302183e7 * Iz3 * t84 - 0.8424409416e4 * T * Iz3 + Iz3 * t145;
      t198 = t5 ^ 2;
      t199 = t198 * t5;
      res__5 = 0.1e1 / t199 / t66 / m1 / LP * (t199 * (t49 * t48 + t7 * (t3 * (t32 * t53 + t41 * t29 - 0.1e1 * t57) + t66 * anc * m1 * LP * t61) + t3 * t76 * anc * t74 + t32 * t81 + t41 * t53 + LP * (t170 + t194)) + t198 * (t49 * t1 * (t32 * t39 + t41 * t81 + t57) + t7 * t1 * (-t76 * anc * t74 + t3 * t48) + t66 * anc * t76 * t3 * t1 * t61) + t49 * (-0.1e1 * t61 * t44 * m1 * Iz3 - 0.1e1 * t61 * m3 * t44 * t64 * m1));
      

      % store on output
      res__h = zeros(5,1);
      res__h(1) = res__1;
      res__h(2) = res__2;
      res__h(3) = res__3;
      res__h(4) = res__4;
      res__h(5) = res__5;

    end
        function res__DhDx = DhDx( self, T, vars__ )

    g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;
      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

      % evaluate function
      t1 = cos(theta1);
      res__1_2 = t1 * LP;
      t2 = cos(theta3);
      res__1_3 = -anc * t2;
      res__2_1 = 1;
      t4 = sin(theta1);
      t5 = t4 * LP;
      res__2_2 = -t5;
      t6 = sin(theta3);
      res__2_3 = anc * t6;
      res__3_2 = -(lambda1 * t1 + lambda2 * t4) * LP;
      res__3_6 = res__2_2;
      res__3_7 = res__1_2;
      t12 = t1 ^ 2;
      res__4_2 = 0.1e1 / t12 * t2 * anc * theta3__dot;
      t19 = 0.1e1 / t1;
      res__4_3 = theta3__dot * t19 * (-t1 * t2 - t6 * t4) * anc;
      res__4_4 = -1;
      t23 = -t6 * t1 + t2 * t4;
      res__4_5 = t19 * t23 * anc;
      t25 = theta3 ^ 2;
      t26 = LP * t25;
      t27 = sin(psi);
      t29 = t27 * LFb * t26;
      t31 = cos(psi);
      t33 = t31 * LFb * t26;
      t36 = t2 ^ 2;
      t42 = t25 * theta3;
      t50 = t12 ^ 2;
      t52 = 0.57404783e5 * t25;
      t53 = g * m3;
      t54 = 0.592621945e4 * t42;
      t55 = 0.182598982e6 * theta3;
      t58 = t31 * LP;
      t60 = 0.474097556e4 * t42;
      t64 = LFb * (-t60 + t25 * 0.490267876799999940e5 + theta3 * (-0.167958569199998077e6) + 0.190914737699991150e6);
      t65 = t27 * LP;
      t68 = LP * lambda2 * anc;
      t72 = 0.1e1 * t53;
      t77 = LP * lambda1 * anc;
      t78 = 0.1e1 * t77;
      t88 = theta3__dot ^ 2;
      t90 = LFb ^ 2;
      t91 = t90 * m3;
      t92 = t91 + Iz3;
      t97 = Iz3 * t88;
      t98 = anc * t97;
      t102 = anc * m3 * t90 * t88;
      t109 = 0.1e1 / LP;
      t110 = 0.1e1 / t92;
      t111 = t110 * t109;
      res__5_2 = 0.1e1 / t50 * t111 * (t50 * (t36 * (0.3e-4 * t29 + 0.8e-4 * t33) + t2 * (t6 * (0.3e-4 * t33 - 0.8e-4 * t29) + LP * (0.3e-6 * t42 + 0.3e-5 * t25))) + t12 * (t36 * (t58 * LFb * (t52 + t53 + 0.1908391204e6 - t54 - t55) + t65 * t64 + t68) + t2 * (t6 * (t58 * t64 + t65 * LFb * (t54 + t55 - 0.1908391204e6 - t52 - t72) - t78) + LP * (0.4266877988e3 * t42 + t25 * (-0.433173389199999929e4) + theta3 * 0.145530395299999982e5 - 0.161877700099999965e5)) + t92 * LP * t6 * t88) + t36 * t4 * (-0.3e1 * t98 - 0.3e1 * t102)) * anc;
      t115 = 0.2316257697e6 * theta3;
      t117 = LFb * (0.645162464e5 * t25 + t53 + 0.274818405e6 - t54 - t115);
      t121 = -t60 + t25 * 0.401374585399999924e5 + theta3 * (-0.110553786099999983e6) + 0.996152466999999888e5;
      t122 = LFb * t121;
      t123 = t65 * t122;
      t127 = 0.6451624645e5 * t25;
      t137 = LP * (0.2133438996e3 * t42 + t25 * (-0.216586694399999988e4) + theta3 * 0.727651976499999819e4 - 0.809388500500000009e4);
      t146 = 0.6400316985e3 * t25 - 0.4331733889e4 * theta3 + 0.5e0 * t97 + 0.5e0 * t90 * m3 * t88 + 0.7276519765e4;
      t149 = 0.5e0 * t53;
      t150 = 0.2963109725e4 * t42;
      t156 = 0.237048778e4 * t42;
      t164 = t12 * t1;
      t206 = 0.1e1 / t164;
      res__5_3 = -0.2000000000e1 * t206 * t111 * anc * (t164 * (t36 * (t58 * t117 + t123 + t68) + t2 * (t6 * (t58 * t122 + t65 * LFb * (t54 + t115 - 0.274818405e6 - t127 - t72) - t78) + t137) + t6 * LP * t146 + t58 * LFb * (-t149 - 0.1793988448e6 + t150 + 0.1403262787e6 * theta3 - 0.3581385488e5 * t25) + t65 * LFb * (-0.415787784999999894e4 + t156 + theta3 * 0.265745014999999948e5 + t25 * (-0.156240646799999959e5)) - 0.5e0 * t68) + t12 * (t36 * t4 * (-t58 * LFb * t121 + t65 * t117 + t77) + t2 * t4 * (t6 * (t58 * LFb * (-t54 + t127 - t115 + 0.274818405e6 + t53) + t123 + t68) - LP * t146) + t4 * (t6 * t137 + t58 * LFb * (-t156 + t25 * 0.245133938600000001e5 + theta3 * (-0.839792845999990241e5) + 0.954573688499955606e5) + t65 * LFb * (-t149 + t150 - 0.2870239154e5 * t25 + 0.91299491e5 * theta3 - 0.954195602e5) - 0.5e0 * t77)) + t2 * t6 * (-0.1e1 * t98 - 0.1e1 * t102));
      res__5_5 = t109 * t206 * theta3__dot * anc * (0.2e1 * t6 * t12 * t5 + 0.2e1 * t2 * t164 * LP - 0.2e1 * t36 * anc);
      t222 = anc ^ 2;
      res__5_6 = t19 * t110 / m1 * (t1 * (-0.1e1 * m1 * t36 * t222 + 0.4e1 * t91 + m1 * t222 + 0.4e1 * Iz3) - 0.1e1 * t2 * t4 * m1 * t6 * t222);
      res__5_7 = 0.1e1 * t19 * t110 * t2 * t23 * t222;
      
      % store on output
      res__DhDx = zeros(5,7);
      res__DhDx(1,2) = res__1_2;
      res__DhDx(1,3) = res__1_3;
      res__DhDx(2,1) = res__2_1;
      res__DhDx(2,2) = res__2_2;
      res__DhDx(2,3) = res__2_3;
      res__DhDx(3,2) = res__3_2;
      res__DhDx(3,6) = res__3_6;
      res__DhDx(3,7) = res__3_7;
      res__DhDx(4,2) = res__4_2;
      res__DhDx(4,3) = res__4_3;
      res__DhDx(4,4) = res__4_4;
      res__DhDx(4,5) = res__4_5;
      res__DhDx(5,2) = res__5_2;
      res__DhDx(5,3) = res__5_3;
      res__DhDx(5,5) = res__5_5;
      res__DhDx(5,6) = res__5_6;
      res__DhDx(5,7) = res__5_7;
    end
    function plot(self, T, vars__)
      s = vars__(1,:);
      theta1 = vars__(2,:);
      theta3 = vars__(3,:);
      s__dot = vars__(4,:);
      theta3__dot = vars__(5,:);
      lambda1 = vars__(6,:);
      lambda2 = vars__(7,:);

      figure(1)
      plot(T,s);
      title("s(t)");

      figure(2)
      plot(T,theta3);
      title("theta3(t)");
 
      figure(3)
      plot(T,theta3__dot);
        title("theta3__dot(t)");

    end
    end
end
    