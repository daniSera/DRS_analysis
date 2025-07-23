classdef DRS_DAE6 < DAC_ODEclass
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
        
        function self = DRS_DAE6(g,m1,m3,psi,anc,LFb,LP,Iz3)
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




function res__f = f( self, t, vars__ )

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

                g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;

      % evaluate function
      res__1 = s__dot;
      t1 = cos(theta3);
      t2 = anc * t1;
      t5 = cos(theta1);
      res__2 = 0.1e1 / t5 / LP * theta3__dot * t2;
      res__3 = theta3__dot;
      res__4 = 0.1e1 / m1 * (lambda1 + act);
      t11 = cos(theta3 + psi);
      t16 = sin(theta3);
      t20 = LFb ^ 2;
      t21 = m3 * t20;
      t22 = t21 + Iz3;
      res__5 = 0.1e1 / t22 * (-g * m3 * t11 * LFb + lambda1 * anc * t16 - lambda2 * t2 + p_Taereo);
      t24 = anc * m1;
      t25 = t5 ^ 2;
      t26 = t25 ^ 2;
      t27 = t26 * t5;
      t28 = sin(theta1);
      t29 = t28 * t27;
      t30 = LP ^ 2;
      t31 = sin(psi);
      t33 = m3 * LFb;
      t34 = t33 * t31 * g;
      t35 = anc * lambda1;
      t36 = t34 + t35;
      t38 = 2 * theta3;
      t39 = sin(t38);
      t41 = t39 * t36 * t30 * t29;
      t43 = cos(psi);
      t45 = t33 * t43 * g;
      t46 = anc * lambda2;
      t47 = t45 + t46;
      t49 = cos(t38);
      t51 = t49 * t47 * t30 * t29;
      t54 = t1 * p_Taereo;
      t55 = theta3__dot ^ 2;
      t60 = t27 * (-t22 * t16 * t55 + t54) * t28 * t30;
      t62 = t47 * LP;
      t63 = t28 * t62;
      t64 = anc * t55;
      t68 = t1 ^ 2;
      t71 = t16 * t36 + p_Taereo;
      t72 = t71 * LP;
      t73 = t1 * t28;
      t76 = LP * t16;
      t78 = t22 * t76 * t28 * t55;
      t79 = t68 * (3 * t22 * t64 + t63) - t73 * t72 + t78;
      t82 = t25 * t1 * t79 * anc;
      t84 = t68 * t1;
      t86 = anc ^ 2;
      t88 = t22 * t86 * t84 * t55;
      t91 = 3 * theta1;
      t92 = cos(t91);
      t94 = sin(t91);
      t98 = LP * t28;
      t105 = 4 * t28 * t22 * t64 - t94 * t28 * t62 + 3 * t5 * t36 * t98 + t25 * t62 - t62;
      t121 = t94 * t36 * t28 / 3 - t25 * t36 / 3 + t5 * t47 * t28 + t34 / 3 + t35 / 3;
      t125 = t22 * t55;
      t127 = p_Taereo * t16;
      t128 = t1 * t125 + t127;
      t131 = t94 * t28 * t26;
      t135 = t47 * t16;
      t137 = t1 * t135 + t68 * t36 - t127 - t34 - t35;
      t139 = t26 ^ 2;
      t141 = lambda1 * t68;
      t146 = t86 * t16 * lambda2 / 2;
      t153 = t26 * t25 * t5;
      t154 = t153 * t1;
      t161 = t86 * t16;
      t163 = t161 * lambda1 * t28 / 2;
      t164 = t36 * LP;
      t168 = t135 + t125 / 8;
      t172 = t34 + 0.7e1 / 0.8e1 * t127 + t35;
      t176 = t26 * t25;
      t181 = t161 * lambda2 * t68;
      t183 = 0.3e1 / 0.4e1 * t54 * t98;
      t184 = 0.3e1 / 0.4e1 * t78;
      t189 = anc * t28;
      t196 = t68 * t28 * anc * t71;
      t200 = t1 * t22 * t55 * LP / 8;
      t202 = p_Taereo * t76 / 8;
      t213 = 0.1e1 / t30;
      t215 = m1 * t86 * t68;
      t216 = m1 * t86;
      t228 = 0.1e1 / (t25 * (t215 - t216 / 2 - t21 / 2 - Iz3 / 2) + m1 * t86 * t5 * t16 * t73 - t215 / 2);
      res__6 = 0.1e1 / t28 / t176 * t228 * t213 * (t92 * (t41 / 8 - t51 / 8 + t60 / 8 + t82 / 8 - 0.3e1 / 0.8e1 * t88) + (t39 * t26 * LP * t105 / 8 - 0.3e1 / 0.8e1 * t49 * t121 * t26 * t30 + t131 * t30 * t128 / 8 + t139 * t137 * t30 - t154 * LP * (t86 * t141 / 2 + t1 * (t63 - t146) - t28 * t72) - t176 * (t86 * lambda2 * t84 * t28 / 2 + t68 * (t163 + t164) + t1 * t168 * LP - t172 * LP) * LP + t27 * (t86 * lambda1 * t84 - t181 + t183 - t184) * LP / 2 + t26 * (t84 * (t45 + 0.3e1 / 0.2e1 * t46) * t189 - t196 + t200 + t202) * LP + 0.3e1 / 0.8e1 * t82 - 0.9e1 / 0.8e1 * t88) * t5) * theta3__dot * t24;
      t245 = t30 * m1;
      t261 = m1 * LP;
      t264 = t21 + t216 + Iz3;
      res__7 = t228 / t153 * t213 * (-t92 * (-t41 + t51 - t60 - t82 + 3 * t88) * m1 / 8 + (t39 * t26 * LP * t105 * m1 / 8 - 0.3e1 / 0.8e1 * t49 * t121 * t26 * t245 + t131 * t30 * m1 * t128 / 8 + t139 * t137 * t245 - t154 * LP * (-t216 * t141 / 2 + t1 * m1 * (t63 + t146) - t28 * t71 * t261 + t264 * lambda1 / 2) - t176 * (-t86 * lambda2 * t84 * t28 * m1 / 2 + t68 * (-t163 + t164) * m1 + t1 * (t28 * t264 * lambda2 / 2 + t168 * t261) - t172 * t261) * LP + t27 * (t181 + t183 - t184) * t261 / 2 + t26 * (t84 * t47 * t189 - t196 + t200 + t202) * t261 + 0.3e1 / 0.8e1 * t25 * t1 * t79 * t24 - 0.9e1 / 0.8e1 * t22 * t86 * t55 * t84 * m1) * t5) * theta3__dot * anc;
      

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
    function res__DfDx = DfDx( self, t, vars__ )

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

                g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;

      % evaluate function
      res__1_1 = s__dot;
      t1 = cos(theta3);
      t2 = anc * t1;
      t5 = cos(theta1);
      res__2_1 = 0.1e1 / t5 / LP * theta3__dot * t2;
      res__3_1 = theta3__dot;
      res__4_1 = 0.1e1 / m1 * (lambda1 + act);
      t11 = cos(theta3 + psi);
      t16 = sin(theta3);
      t20 = LFb ^ 2;
      t21 = m3 * t20;
      t22 = t21 + Iz3;
      res__5_1 = 0.1e1 / t22 * (-g * m3 * t11 * LFb + lambda1 * anc * t16 - lambda2 * t2 + p_Taereo);
      t24 = anc * m1;
      t25 = t5 ^ 2;
      t26 = t25 ^ 2;
      t27 = t26 * t5;
      t28 = sin(theta1);
      t29 = t28 * t27;
      t30 = LP ^ 2;
      t31 = sin(psi);
      t33 = m3 * LFb;
      t34 = t33 * t31 * g;
      t35 = anc * lambda1;
      t36 = t34 + t35;
      t38 = 2 * theta3;
      t39 = sin(t38);
      t41 = t39 * t36 * t30 * t29;
      t43 = cos(psi);
      t45 = t33 * t43 * g;
      t46 = anc * lambda2;
      t47 = t45 + t46;
      t49 = cos(t38);
      t51 = t49 * t47 * t30 * t29;
      t54 = t1 * p_Taereo;
      t55 = theta3__dot ^ 2;
      t60 = t27 * (-t22 * t16 * t55 + t54) * t28 * t30;
      t62 = t47 * LP;
      t63 = t28 * t62;
      t64 = anc * t55;
      t68 = t1 ^ 2;
      t71 = t16 * t36 + p_Taereo;
      t72 = t71 * LP;
      t73 = t1 * t28;
      t76 = LP * t16;
      t78 = t22 * t76 * t28 * t55;
      t79 = t68 * (3 * t22 * t64 + t63) - t73 * t72 + t78;
      t82 = t25 * t1 * t79 * anc;
      t84 = t68 * t1;
      t86 = anc ^ 2;
      t88 = t22 * t86 * t84 * t55;
      t91 = 3 * theta1;
      t92 = cos(t91);
      t94 = sin(t91);
      t98 = LP * t28;
      t105 = 4 * t28 * t22 * t64 - t94 * t28 * t62 + 3 * t5 * t36 * t98 + t25 * t62 - t62;
      t121 = t94 * t36 * t28 / 3 - t25 * t36 / 3 + t5 * t47 * t28 + t34 / 3 + t35 / 3;
      t125 = t22 * t55;
      t127 = p_Taereo * t16;
      t128 = t1 * t125 + t127;
      t131 = t94 * t28 * t26;
      t135 = t47 * t16;
      t137 = t1 * t135 + t68 * t36 - t127 - t34 - t35;
      t139 = t26 ^ 2;
      t141 = lambda1 * t68;
      t146 = t86 * t16 * lambda2 / 2;
      t153 = t26 * t25 * t5;
      t154 = t153 * t1;
      t161 = t86 * t16;
      t163 = t161 * lambda1 * t28 / 2;
      t164 = t36 * LP;
      t168 = t135 + t125 / 8;
      t172 = t34 + 0.7e1 / 0.8e1 * t127 + t35;
      t176 = t26 * t25;
      t181 = t161 * lambda2 * t68;
      t183 = 0.3e1 / 0.4e1 * t54 * t98;
      t184 = 0.3e1 / 0.4e1 * t78;
      t189 = anc * t28;
      t196 = t68 * t28 * anc * t71;
      t200 = t1 * t22 * t55 * LP / 8;
      t202 = p_Taereo * t76 / 8;
      t213 = 0.1e1 / t30;
      t215 = m1 * t86 * t68;
      t216 = m1 * t86;
      t228 = 0.1e1 / (t25 * (t215 - t216 / 2 - t21 / 2 - Iz3 / 2) + m1 * t86 * t5 * t16 * t73 - t215 / 2);
      res__6_1 = 0.1e1 / t28 / t176 * t228 * t213 * (t92 * (t41 / 8 - t51 / 8 + t60 / 8 + t82 / 8 - 0.3e1 / 0.8e1 * t88) + (t39 * t26 * LP * t105 / 8 - 0.3e1 / 0.8e1 * t49 * t121 * t26 * t30 + t131 * t30 * t128 / 8 + t139 * t137 * t30 - t154 * LP * (t86 * t141 / 2 + t1 * (t63 - t146) - t28 * t72) - t176 * (t86 * lambda2 * t84 * t28 / 2 + t68 * (t163 + t164) + t1 * t168 * LP - t172 * LP) * LP + t27 * (t86 * lambda1 * t84 - t181 + t183 - t184) * LP / 2 + t26 * (t84 * (t45 + 0.3e1 / 0.2e1 * t46) * t189 - t196 + t200 + t202) * LP + 0.3e1 / 0.8e1 * t82 - 0.9e1 / 0.8e1 * t88) * t5) * theta3__dot * t24;
      t245 = t30 * m1;
      t261 = m1 * LP;
      t264 = t21 + t216 + Iz3;
      res__7_1 = t228 / t153 * t213 * (-t92 * (-t41 + t51 - t60 - t82 + 3 * t88) * m1 / 8 + (t39 * t26 * LP * t105 * m1 / 8 - 0.3e1 / 0.8e1 * t49 * t121 * t26 * t245 + t131 * t30 * m1 * t128 / 8 + t139 * t137 * t245 - t154 * LP * (-t216 * t141 / 2 + t1 * m1 * (t63 + t146) - t28 * t71 * t261 + t264 * lambda1 / 2) - t176 * (-t86 * lambda2 * t84 * t28 * m1 / 2 + t68 * (-t163 + t164) * m1 + t1 * (t28 * t264 * lambda2 / 2 + t168 * t261) - t172 * t261) * LP + t27 * (t181 + t183 - t184) * t261 / 2 + t26 * (t84 * t47 * t189 - t196 + t200 + t202) * t261 + 0.3e1 / 0.8e1 * t25 * t1 * t79 * t24 - 0.9e1 / 0.8e1 * t22 * t86 * t55 * t84 * m1) * t5) * theta3__dot * anc;
      
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
    function res__h = h( self, t, vars__ )

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

                g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;

      % evaluate function
      t1 = sin(theta1);
      t2 = t1 * LP;
      t3 = sin(theta3);
      t4 = anc * t3;
      res__1 = t2 - t4 + y1 - y4;
      t5 = cos(theta1);
      t7 = cos(theta3);
      res__2 = t5 * LP - anc * t7 + L1 + s + x1 - x4;
      res__3 = (-lambda1 * t1 + lambda2 * t5) * LP;
      res__4 = -s__dot + theta3__dot / t5 * (t1 * t7 - t5 * t3) * anc;
      t19 = anc * m1;
      t20 = sin(psi);
      t21 = t20 * g;
      t22 = m3 * LFb;
      t25 = anc * lambda1 + t22 * t21;
      t26 = t7 ^ 2;
      t29 = cos(psi);
      t33 = t22 * t29 * g + anc * lambda2;
      t35 = theta3__dot ^ 2;
      t36 = LFb ^ 2;
      t38 = m3 * t36 + Iz3;
      t48 = anc ^ 2;
      t55 = t5 ^ 2;
      t56 = t55 * t5;
      res__5 = 0.1e1 / t38 / m1 / LP / t56 * (-t56 * (t26 * t25 * t19 + t7 * (t33 * t3 - t38 * t35) * t19 - t4 * p_Taereo * m1 - anc * t22 * m1 * t21 - m1 * lambda1 * t48 - (lambda1 + act) * t38) * LP + t55 * t2 * anc * (t26 * t33 + t7 * (-t3 * t25 - p_Taereo) + t38 * t3 * t35) * m1 - t38 * t48 * m1 * t26 * t35);
      

      % store on output
      res__h = zeros(5,1);
      res__h(1) = res__1;
      res__h(2) = res__2;
      res__h(3) = res__3;
      res__h(4) = res__4;
      res__h(5) = res__5;

    end
    function res__DhDx = DhDx( self, t, vars__ )

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta3 = vars__(3);
      s__dot = vars__(4);
      theta3__dot = vars__(5);
      lambda1 = vars__(6);
      lambda2 = vars__(7);

                g   = self.g  ;
    m1  = self.m1 ;
    m3  = self.m3 ;
    psi = self.psi;
    anc = self.anc;
    LFb = self.LFb;
    LP  = self.LP ;
    Iz3 = self.Iz3;
    
      % evaluate function
      t1 = cos(theta1);
      res__1_2 = t1 * LP;
      t2 = cos(theta3);
      t3 = anc * t2;
      res__1_3 = -t3;
      res__2_1 = 1;
      t4 = sin(theta1);
      t5 = t4 * LP;
      res__2_2 = -t5;
      t6 = sin(theta3);
      res__2_3 = anc * t6;
      res__3_2 = -(lambda1 * t1 + t4 * lambda2) * LP;
      res__3_6 = res__2_2;
      res__3_7 = res__1_2;
      t11 = t1 ^ 2;
      res__4_2 = 0.1e1 / t11 * theta3__dot * t3;
      t18 = 0.1e1 / t1;
      res__4_3 = theta3__dot * t18 * (-t1 * t2 - t4 * t6) * anc;
      res__4_4 = -1;
      t20 = t4 * t2;
      t22 = -t1 * t6 + t20;
      res__4_5 = t18 * t22 * anc;
      t24 = cos(psi);
      t26 = m3 * LFb;
      t27 = t26 * t24 * g;
      t28 = anc * lambda2;
      t29 = t27 + t28;
      t30 = t2 ^ 2;
      t31 = t30 * t29;
      t32 = sin(psi);
      t34 = t26 * t32 * g;
      t35 = anc * lambda1;
      t36 = -t34 - t35;
      t37 = t6 * t36;
      t40 = theta3__dot ^ 2;
      t42 = LFb ^ 2;
      t43 = m3 * t42;
      t44 = t43 + Iz3;
      t45 = t44 * t6 * t40;
      t56 = t11 ^ 2;
      t58 = 0.1e1 / LP;
      t60 = 0.1e1 / t44;
      res__5_2 = t60 * t58 / t56 * anc * (t11 * LP * (t31 + t2 * (t37 - p_Taereo) + t45) - 3 * t44 * anc * t4 * t30 * t40);
      t70 = t11 * t1;
      t94 = 0.1e1 / t70;
      res__5_3 = -2 * t94 * t60 * t58 * anc * (t70 * (t31 + t2 * (t37 - p_Taereo / 2) + t45 / 2 - t27 / 2 - t28 / 2) * LP + t11 * t4 * LP * (-t30 * t36 + t2 * (t29 * t6 - t44 * t40 / 2) - t34 / 2 - t35 / 2 - p_Taereo * t6 / 2) - t44 * anc * t6 * t2 * t40);
      res__5_5 = 2 * t94 * (t70 * t2 * LP + t6 * t11 * t5 - t30 * anc) * theta3__dot * anc * t58;
      t108 = anc ^ 2;
      t120 = t60 * t18;
      res__5_6 = t120 / m1 * (t1 * (-m1 * t108 * t30 + m1 * t108 + Iz3 + t43) - m1 * t108 * t6 * t20);
      res__5_7 = t2 * t120 * t22 * t108;
      
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
