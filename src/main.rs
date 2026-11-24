use plotters::prelude::*;

mod tabela_ar;
use std::sync::Arc;
use std::time::Instant;
use tabela_ar::TabelaAr;

fn main() {
    println!("Hello, world!");

    //tabela de prop fisicas do ar:::
    let tabela_ar = TabelaAr::nova(101_325.0).expect("Erro ao criar tabela de propriedades do ar");

    // Exemplo: ar a 350 °C
    let props = tabela_ar.consultar(350.0);

    println!("Densidade: {}", props.rho);
    println!("Calor específico: {}", props.cp);
    println!("Condutividade: {}", props.k);
    println!("Viscosidade: {}", props.mu);

    //definindo a malha de pontos
    let tamanho: f64 = 1.0;
    const N: usize = 1000;
    const ITERACOES_TEMPO: usize = 1_000_000;
    let n: f64 = N as f64;
    let dx = tamanho / n;
    let dt = 0.003;
    let tabela_ar =
        Arc::new(TabelaAr::nova(101_325.0).expect("Erro ao inicializar propriedades do ar"));

    //condições de contorno
    let ta = 0.0;
    let tb = 350.;

    //meios
    //TODO: que tal criar malhas com valores variados de temperatura pra estressar os métodos??
    //cria os materiais
    let mut meio1 = Meio::criar_material(N / 2, 350.0, Material::Ferro);
    let mut meio2 = Meio::criar_material(N / 2, 350.0, Material::Ar(Arc::clone(&tabela_ar)));

    //define o valor da velocidade em cada nó
    let pulando_solido = meio1.temperaturas.len() as f64 * dx;
    Meio::criar_vetor_advectivo(
        &mut meio1.velocidade_advectiva,
        |x: f64| meio1.material.get_u(pulando_solido + x),
        dx,
    );
    Meio::criar_vetor_advectivo(
        &mut meio2.velocidade_advectiva,
        |x: f64| meio2.material.get_u(pulando_solido + x),
        dx,
    );
    //let mut meio2 = Meio::criar_material(N / 2, 350.0, Material::Cobre);

    meio1.temperaturas[0] = ta;
    //meio2.temperaturas[N / 2 - 1] = tb;
    let ultimo = meio2.temperaturas.len() - 1;
    meio2.temperaturas[ultimo] = tb;

    //metodo
    let metodinho = Metodo {
        dimensao: tamanho,
        nos: N,
        dx: dx,
        dt: dt,
    };

    //declarando o sistema
    let mut sistema_impl = Sistema {
        interface: meio1.temperaturas.len(),
        meio1: meio1.clone(),
        meio2: meio2.clone(),
        metodo: metodinho.clone(),
    };

    let mut sistema_expl = Sistema {
        interface: meio1.temperaturas.len(),
        meio1,
        meio2,
        metodo: metodinho,
    };

    //teste de fluxo
    //let fluxo = interface_flux(&sistema, 0);
    //println!("Fluxo entre os nós 0 e 1: {:?}", fluxo);

    //teste que eu tenho que arrumar depois
    //let mut novo_meio1 = sistema.meio1.temperaturas.clone();
    //let mut novo_meio2 = sistema.meio2.temperaturas.clone();

    //sistema.iteracao(&mut novo_meio1, &mut novo_meio2);

    //println!("{:?}", &novo_meio1[495..500]);
    //println!("{:?}", &novo_meio2[495..500]);

    //simulação implicita
    //sistema_expl.calc_explicito(ITERACOES_TEMPO);
    println!("metodo explicito:");
    benchmark(ITERACOES_TEMPO, || {
        sistema_expl.calc_explicito(ITERACOES_TEMPO);
    });

    let tmax = sistema_expl
        .meio1
        .temperaturas
        .iter()
        .chain(sistema_expl.meio2.temperaturas.iter())
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    println!("Temperatura máxima: {tmax:.3} °C");

    sistema_expl.plotar("explicito.png", "Método explícito");

    //parte do implicito
    let mut tridiag =
        SistemaTridiagonal::criar(N.checked_sub(2).expect("N tem que ser maior que 2!!"));
    let coeficientes = sistema_impl.coef_implicit_adv(6);
    println!("Coeficientes do nó 6: {:?}", coeficientes);
    //sistema_impl.met_implicit(&mut tridiag);
    let mut sis_impl2 = sistema_impl.clone();
    let mut sis_impl3 = sistema_impl.clone();

    println!("metodo implicito twisted:");
    benchmark(ITERACOES_TEMPO, || {
        sistema_impl.solver_implicit_claude(&mut tridiag, ITERACOES_TEMPO);
    });

    println!("metodo implicito homemade:");
    benchmark(ITERACOES_TEMPO, || {
        sis_impl2.solver_implicit(&mut tridiag, ITERACOES_TEMPO);
    });

    println!("metodo implicito classico:");
    benchmark(ITERACOES_TEMPO, || {
        sis_impl3.solver_implicit_classico(&mut tridiag, ITERACOES_TEMPO);
    });
    //loop completo
    // let mut tridiag =
    //     SistemaTridiagonal::criar(N.checked_sub(2).expect("N tem que ser maior que 2!!"));
    //
    // println!("len b = {}", tridiag.b.len());
    // println!(
    //     "destino meio1 = {}",
    //     sistema_impl.meio1.temperaturas[1..].len()
    // );
    //
    // println!("{:?}", &tridiag.diagonal_principal[..5]);
    // println!("{:?}", &tridiag.b[..5]);
    //
    // sistema_impl.solver_implicit(&mut tridiag, ITERACOES_TEMPO);
    //
    // //vetor completo
    // let mut temperaturas = Vec::with_capacity(N);
    //
    // temperaturas.extend_from_slice(&sistema_impl.meio1.temperaturas);
    // temperaturas.extend_from_slice(&sistema_impl.meio2.temperaturas);
    //
    comparacao(&sistema_impl, &sistema_expl);
    maior_diferenca(&sistema_impl, &sistema_expl);
    //

    println!("diferença entre os mesmos métodos:");
    comparacao(&sistema_impl, &sis_impl2);
    maior_diferenca(&sistema_impl, &sis_impl2);

    println!("diferença entre thomas classico e o meu homemade:");
    comparacao(&sis_impl3, &sis_impl2);
    maior_diferenca(&sis_impl3, &sis_impl2);

    let tmax = sistema_impl
        .meio1
        .temperaturas
        .iter()
        .chain(sistema_impl.meio2.temperaturas.iter())
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    println!("Temperatura máxima implícita: {tmax:.3} °C");

    println!(
        "u no início do ar: {}",
        sistema_expl.meio2.velocidade_advectiva[0]
    );
    println!(
        "u no meio do ar: {}",
        sistema_expl.meio2.velocidade_advectiva[250]
    );

    //plot
    sistema_impl.plotar("implicito.png", "Método implícito");
}

#[derive(Clone)]
struct Meio {
    temperaturas: Vec<f64>,
    velocidade_advectiva: Vec<f64>,
    material: Material,
}

#[derive(Clone)]
struct Metodo {
    dimensao: f64,
    nos: usize,
    dx: f64,
    dt: f64,
}

#[derive(Clone)]
struct Sistema {
    meio1: Meio,
    meio2: Meio,
    metodo: Metodo,
    interface: usize,
}

struct SistemaTridiagonal {
    diagonal_inferior: Vec<f64>,
    diagonal_principal: Vec<f64>,
    diagonal_superior: Vec<f64>,
    b: Vec<f64>,
}

impl SistemaTridiagonal {
    fn criar(tamanho: usize) -> Self {
        Self {
            diagonal_inferior: vec![0.0; tamanho],
            diagonal_principal: vec![0.0; tamanho],
            diagonal_superior: vec![0.0; tamanho],
            b: vec![0.0; tamanho],
        }
    }
}

#[derive(Clone)]
enum Material {
    Ferro,
    Cobre,
    Ar(Arc<TabelaAr>),
}

impl Material {
    fn get_k(&self, t: f64) -> f64 {
        match self {
            Material::Ferro => 10.0 + 0.5 * t,
            Material::Cobre => 15.0 + 0.3 * t,
            Material::Ar(tabela) => tabela.consultar(t).k,
        }
    }

    fn get_rho(&self, t: f64) -> f64 {
        match self {
            Material::Ferro => 7874.0,
            Material::Cobre => 8960.0,
            Material::Ar(tabela) => tabela.consultar(t).rho,
        }
    }

    fn get_cp(&self, t: f64) -> f64 {
        match self {
            Material::Ferro => 449.0,
            Material::Cobre => 385.0,
            Material::Ar(tabela) => tabela.consultar(t).cp,
        }
    }
    fn get_u(&self, x: f64) -> f64 {
        match self {
            Material::Ar(_) => {
                let x_interface = 0.5;
                let l = 0.02;
                let u = -0.1;

                let d = (x - x_interface).max(0.0);
                u * (1.0 - (-d / l).exp())
            }
            _ => 0.0,
        }
    }
}

impl Meio {
    fn criar_material(tamanho: usize, valor: f64, material: Material) -> Self {
        Self {
            temperaturas: vec![valor; tamanho],
            material: material,
            velocidade_advectiva: vec![0.0; tamanho],
        }
    }
    fn criar_vetor_advectivo<F>(velocidades: &mut [f64], f: F, dx: f64)
    where
        F: Fn(f64) -> f64,
    {
        // j: índice do nó, começando em zero.
        // x_no = (j + 0.5) * dx: posição do centro do nó.
        // x_face = (j + 1.0) * dx: posição da face direita desse nó.
        for j in 0..velocidades.len() {
            // 1. Calcula x_j
            // 2. Calcula u_j chamando f
            velocidades[j] = f((j as f64 + 0.5) * dx);
            // 3. Armazena u_j no vetor}
        }
    }
}

impl Sistema {
    fn indice(&self, posicao: usize) -> (&Meio, usize) {
        let q1 = self.meio1.temperaturas.len();

        if posicao < q1 {
            return (&self.meio1, posicao);
        } else {
            return (&self.meio2, posicao - q1);
        }
    }

    fn iteracao(&self, buffer1: &mut Vec<f64>, buffer2: &mut Vec<f64>) {
        //trocar os crone por memory swap
        //não é tão simples assim
        //a pira aqui seria passar dois vetores já criados e modificar eles ou então modificar
        //direto as temperaturas dos meios que são passadas no self.

        //std::mem::swap(&mut t_old, &mut t_new); //colouqei aqui só pra dar ctrl c rápidor

        let n = self.metodo.nos;

        let q1 = self.meio1.temperaturas.len();

        let (mut q_esquerda, _) = interface_flux(self, 0);

        for j in 1..n - 1 {
            let (mei, pos) = self.indice(j);
            //fazer um esquema de reaproveitar o fluxo da direita na esquerda do próximo pra não
            //precisar recalcular
            let (q_direita, _) = interface_flux(self, j);

            let atualtemp = mei.temperaturas[pos];

            let rho = mei.material.get_rho(atualtemp);
            let cp = mei.material.get_cp(atualtemp);
            let termo_sei_la_o_nome = self.metodo.dt / (rho * cp * self.metodo.dx);

            let termo_advectivo = self.adveccao_no(j);

            let novatemp = atualtemp + termo_sei_la_o_nome * (q_esquerda - q_direita)
                - termo_advectivo * self.metodo.dt;
            q_esquerda = q_direita;

            if j < q1 {
                buffer1[pos] = novatemp
            } else {
                buffer2[pos] = novatemp
            }
        }
    }

    //nao usamos isso aqui pelo menos por enquanto
    //os indices das faces precisam ser tratados e se adequar a lógica do sistema.
    fn adveccao_face(&self, face: usize) -> f64 {
        let inicio_ar = self.meio1.temperaturas.len();

        // Sem advecção no sólido nem na interface sólido-fluido
        if face < inicio_ar {
            return 0.0;
        }
        let (meio, pos_local) = self.indice(face);

        // let j = face - inicio_ar;
        //
        // let x = (face + 1) as f64 * self.metodo.dx;
        //
        // let u = self.meio2.material.get_u(x);
        //
        //u da face tem que ser calculado pela media dos nós pra fazer sentido ter a malha de u!
        let u =
            (meio.velocidade_advectiva[pos_local] + meio.velocidade_advectiva[pos_local + 1]) * 0.5;

        //caso muito específico mas que pode economizar processamento
        if u == 0.0 {
            return 0.0;
        }

        // Propriedades congeladas no instante anterior
        let temp_upwind = if u >= 0.0 {
            self.meio2.temperaturas[pos_local]
        } else {
            self.meio2.temperaturas[pos_local + 1]
        };

        let cp = self.meio2.material.get_cp(temp_upwind);
        let rho_ar = meio.material.get_rho(temp_upwind);

        rho_ar * cp * u
    }

    // u dt/dx bacaninha mais fácil de usar
    fn adveccao_no(&self, posicao_global: usize) -> f64 {
        // Não existe advecção no sólido
        if posicao_global < self.interface {
            return 0.0;
        }
        let (blebleble, posicao) = self.indice(posicao_global);
        let u = blebleble.velocidade_advectiva[posicao];
        if u == 0.0 || (posicao == 0 && u > 0.0) {
            return 0.0;
        }
        let dx = self.metodo.dx;
        if u < 0.0 && posicao + 1 >= blebleble.temperaturas.len() {
            panic!(
                "Fora do vetor! global={}, local={}, len={}, u={}",
                posicao_global,
                posicao,
                blebleble.temperaturas.len(),
                u
            );
        }
        let dt = if u > 0.0 {
            blebleble.temperaturas[posicao] - blebleble.temperaturas[posicao - 1]
        } else {
            blebleble.temperaturas[posicao + 1] - blebleble.temperaturas[posicao]
        };
        //let sinal = if u >= 0.0 { 1.0 } else { -1.0 };

        u * dt / dx
    }

    fn coef_implicit_adv(&self, indice_global: usize) -> (f64, f64, f64, f64) {
        let (meio, indice_local) = self.indice(indice_global);

        let temp_atual = meio.temperaturas[indice_local];
        let rho_ar = meio.material.get_rho(temp_atual);

        let rho = if indice_global >= self.meio1.temperaturas.len() {
            rho_ar
        } else {
            meio.material.get_rho(temp_atual)
        };

        let cp = meio.material.get_cp(temp_atual);
        let capacidade = (rho * cp * self.metodo.dx) / self.metodo.dt;

        let (_, ge) = interface_flux(self, indice_global - 1);
        let (_, gd) = interface_flux(self, indice_global);

        // Coeficientes advectivos nas duas interfaces
        let fe = self.adveccao_face(indice_global - 1);
        let fd = self.adveccao_face(indice_global);

        let aw = -ge - fe.max(0.0);
        let ap = capacidade + ge + gd + fd.max(0.0) + (-fe).max(0.0) - (fd - fe);
        let ae = -gd - (-fd).max(0.0);
        let b = capacidade * temp_atual;

        (aw, ap, ae, b)
    }

    //essa função só vai retornar os coeficientes pra usar em cada posição da matriz
    fn coef_implicit(&self, indice_global: usize) -> (f64, f64, f64, f64) {
        let (meio, indice_local) = self.indice(indice_global);

        let temp_atual = meio.temperaturas[indice_local];

        let rho = meio.material.get_rho(temp_atual);
        let cp = meio.material.get_cp(temp_atual);
        let capacidade = (rho * cp * self.metodo.dx) / self.metodo.dt;

        let (_, ge) = interface_flux(self, indice_global - 1);
        let (_, gd) = interface_flux(self, indice_global);

        //aW, aP, aE, b
        (-ge, capacidade + ge + gd, -gd, capacidade * temp_atual)
    }

    fn calc_explicito(&mut self, num_iteracoes: usize) {
        let mut novas_temperaturas_meio1 = self.meio1.temperaturas.clone();
        let mut novas_temperaturas_meio2 = self.meio2.temperaturas.clone();

        for _ in 0..num_iteracoes {
            self.iteracao(&mut novas_temperaturas_meio1, &mut novas_temperaturas_meio2);

            std::mem::swap(&mut self.meio1.temperaturas, &mut novas_temperaturas_meio1);
            std::mem::swap(&mut self.meio2.temperaturas, &mut novas_temperaturas_meio2);

            //colocar um if nesse println pra não fazer 1 MELHAO de prints!!
            // let passo = i + 1;
            // if passo % 100000 == 0 {
            //     println!(
            //         "iteração {:?}, começo: {:?} fim: {:?}",
            //         i,
            //         &self.meio1.temperaturas[0..5],
            //         &self.meio2.temperaturas[495..500],
            //     );
            //     //print da interface
            //     println!("Fe: {:?}", &self.meio1.temperaturas[495..500]);
            //     println!("Cu: {:?}", &self.meio2.temperaturas[0..5]);
            // }
        }
    }

    //função que monta a matriz ou o sistema de equações pra ser resolvido
    fn met_implicit(&self, sistema: &mut SistemaTridiagonal) {
        let tamanho = self
            .metodo
            .nos
            .checked_sub(2)
            .expect("o tamanho da malha deve ser maior que 2!!!!");

        for i in 0..tamanho {
            (
                sistema.diagonal_inferior[i],
                sistema.diagonal_principal[i],
                sistema.diagonal_superior[i],
                sistema.b[i],
            ) = self.coef_implicit_adv(i + 1);
        }

        let ultima_linha = tamanho - 1;
        let ponta_direita = self.meio2.temperaturas.len() - 1;
        //tratamento primeira linha
        sistema.b[0] = sistema.b[0] - sistema.diagonal_inferior[0] * self.meio1.temperaturas[0];
        sistema.diagonal_inferior[0] = 0.0;
        //tratamento ultima linha
        sistema.b[ultima_linha] = sistema.b[ultima_linha]
            - sistema.diagonal_superior[ultima_linha] * self.meio2.temperaturas[ponta_direita];
        sistema.diagonal_superior[ultima_linha] = 0.0;
    }

    //algoritmo de thomas
    fn solver_implicit(&mut self, mut sistema: &mut SistemaTridiagonal, num_iteracoes: usize) {
        let n1 = self.meio1.temperaturas.len();
        let n2 = self.meio2.temperaturas.len();

        for _ in 0..num_iteracoes {
            self.met_implicit(&mut sistema);
            for i in 0..sistema.diagonal_principal.len() - 1 {
                let c1 = sistema.diagonal_inferior[i + 1] / sistema.diagonal_principal[i];
                sistema.diagonal_inferior[i] = 0.0;
                sistema.diagonal_principal[i] *= -c1;
                sistema.diagonal_superior[i] *= -c1;
                sistema.b[i] *= -c1;

                sistema.diagonal_inferior[i + 1] += sistema.diagonal_principal[i];
                sistema.diagonal_principal[i + 1] += sistema.diagonal_superior[i];
                //diagonal_superior[i + 1] = diagonal_superior[i +1];
                sistema.b[i + 1] += sistema.b[i];
            }
            let ultimo_termo = sistema.b.len() - 1;
            sistema.b[ultimo_termo] =
                sistema.b[ultimo_termo] / sistema.diagonal_principal[ultimo_termo];
            sistema.diagonal_principal[ultimo_termo] = 1.0;
            for n in (1..sistema.diagonal_principal.len()).rev() {
                //agora o negócio é pegar o bn e dividir por pn
                //salva o tn
                //faz bn-1 = bn-1 - sn-1 * tn
                //divide bn por pn-1 e salva no tn-1
                sistema.b[n - 1] = (sistema.b[n - 1]
                    - sistema.diagonal_superior[n - 1] * sistema.b[n])
                    / sistema.diagonal_principal[n - 1];
                sistema.diagonal_principal[n - 1] = 1.0;

                //repete até o indice zero
            }
            //atualiza os dados das temperaturas com os b obtidos
            self.meio1.temperaturas[1..].copy_from_slice(&sistema.b[..n1 - 1]);
            self.meio2.temperaturas[..n2 - 1].copy_from_slice(&sistema.b[n1 - 1..]);
        }
    }

    //feito por ia a fim de otimizar o tempo
    fn solver_implicit_classico(
        &mut self,
        _sistema: &mut SistemaTridiagonal,
        num_iteracoes: usize,
    ) {
        let n1 = self.meio1.temperaturas.len();
        let n2 = self.meio2.temperaturas.len();
        let n = n1 + n2;
        assert!(n1 >= 1 && n2 >= 1 && n >= 4);

        let dx = self.metodo.dx;
        let cap_fator = dx / self.metodo.dt;
        let dois_sobre_dx = 2.0 / dx;

        // Temperaturas em um vetor contínuo
        let mut t = Vec::with_capacity(n);
        t.extend_from_slice(&self.meio1.temperaturas);
        t.extend_from_slice(&self.meio2.temperaturas);

        // Buffers alocados uma única vez
        let mut k = vec![0.0; n];
        let mut rc = vec![0.0; n];
        let mut g = vec![0.0; n - 1];
        let mut fl = vec![0.0; n - 1];

        let mut aw = vec![0.0; n];
        let mut ap = vec![0.0; n];
        let mut ae = vec![0.0; n];
        let mut bb = vec![0.0; n];

        let mut c_l = vec![0.0; n];
        let mut d_l = vec![0.0; n];

        // Velocidades nas faces
        let mut u_face = vec![0.0; n - 1];
        let v = &self.meio2.velocidade_advectiva;

        for f in n1..n - 1 {
            u_face[f] = (v[f - n1] + v[f - n1 + 1]) * 0.5;
        }

        for _ in 0..num_iteracoes {
            // 1. Propriedades termofísicas
            preencher_props(&self.meio1.material, &t[..n1], &mut k[..n1], &mut rc[..n1]);

            preencher_props(&self.meio2.material, &t[n1..], &mut k[n1..], &mut rc[n1..]);

            // 2. Condutâncias e fluxos advectivos
            for (f, (g, fl)) in g.iter_mut().zip(fl.iter_mut()).enumerate() {
                let (ka, kb) = (k[f], k[f + 1]);

                *g = dois_sobre_dx * ka * kb / (ka + kb);

                let u = u_face[f];
                *fl = if u >= 0.0 { rc[f] } else { rc[f + 1] } * u;
            }

            // 3. Montagem do sistema tridiagonal
            for i in 1..n - 1 {
                let (ge, gd) = (g[i - 1], g[i]);
                let (fe, fd) = (fl[i - 1], fl[i]);

                let cap = rc[i] * cap_fator;

                aw[i] = -ge - fe.max(0.0);

                ap[i] = cap + ge + gd + fd.max(0.0) + (-fe).max(0.0) - (fd - fe);

                ae[i] = -gd - (-fd).max(0.0);

                bb[i] = cap * t[i];
            }

            // 4. THOMAS CLÁSSICO: eliminação direta
            let mut cc = 0.0;
            let mut dd = t[0];

            for i in 1..n - 1 {
                let inv = 1.0 / (ap[i] - aw[i] * cc);

                cc = ae[i] * inv;
                dd = (bb[i] - aw[i] * dd) * inv;

                c_l[i] = cc;
                d_l[i] = dd;
            }

            // 5. THOMAS CLÁSSICO: substituição regressiva
            // t[n-1] permanece fixo (condição de contorno)
            for i in (1..n - 1).rev() {
                t[i] = d_l[i] - c_l[i] * t[i + 1];
            }
        }

        // Atualiza temperaturas
        self.meio1.temperaturas.copy_from_slice(&t[..n1]);
        self.meio2.temperaturas.copy_from_slice(&t[n1..]);
    }

    // ↓ substitui o solver_implicit antigo
    /// Implícito com advecção: montagem + Thomas fundidos, sem alocação no loop.
    /// `_tridiag` não é mais usado (fica só pra não mexer no main).
    fn solver_implicit_claude(&mut self, _tridiag: &mut SistemaTridiagonal, num_iteracoes: usize) {
        let n1 = self.meio1.temperaturas.len();
        let n2 = self.meio2.temperaturas.len();
        let n = n1 + n2;
        assert!(n1 >= 1 && n2 >= 1 && n >= 4);

        // linhas internas 1..=n-2, divididas em metade de cima (1..=m) e de baixo (m+1..=n-2)
        let r = n - 2;
        let h = r / 2;
        let m = r - h;

        let dx = self.metodo.dx;
        let cap_fator = dx / self.metodo.dt;
        let dois_sobre_dx = 2.0 / dx;

        // estado contínuo: as temperaturas ficam num vetor só durante toda a simulação
        let mut t = Vec::with_capacity(n);
        t.extend_from_slice(&self.meio1.temperaturas);
        t.extend_from_slice(&self.meio2.temperaturas);

        // buffers, alocados uma vez
        let mut k = vec![0.0; n]; // condutividade por nó
        let mut rc = vec![0.0; n]; // rho*cp por nó
        let mut g = vec![0.0; n - 1]; // condutância da face f (entre nós f e f+1)
        let mut fl = vec![0.0; n - 1]; // rho*cp*u upwind na face f
        let mut aw = vec![0.0; n];
        let mut ap = vec![0.0; n];
        let mut ae = vec![0.0; n];
        let mut bb = vec![0.0; n];
        let mut c_l = vec![0.0; n]; // coeficientes eliminados (c' em cima, e' embaixo)
        let mut d_l = vec![0.0; n]; // lado direito eliminado (d' em cima, f' embaixo)

        // u nas faces não depende do tempo: calcula uma vez (zero no sólido e na interface)
        let mut u_face = vec![0.0; n - 1];
        let v = &self.meio2.velocidade_advectiva;
        for f in n1..n - 1 {
            u_face[f] = (v[f - n1] + v[f - n1 + 1]) * 0.5;
        }

        for _ in 0..num_iteracoes {
            // 1) propriedades por nó (cada nó avaliado uma vez só)
            preencher_props(&self.meio1.material, &t[..n1], &mut k[..n1], &mut rc[..n1]);
            preencher_props(&self.meio2.material, &t[n1..], &mut k[n1..], &mut rc[n1..]);

            // 2) faces: cada face calculada uma vez e usada pelos dois nós vizinhos.
            //    g = 1/(dx/2k1 + dx/2k2) = 2 k1 k2 / (dx (k1+k2)) → uma divisão em vez de três
            for (f, (g, fl)) in g.iter_mut().zip(fl.iter_mut()).enumerate() {
                let (ka, kb) = (k[f], k[f + 1]);
                *g = dois_sobre_dx * ka * kb / (ka + kb);
                let u = u_face[f];
                *fl = if u >= 0.0 { rc[f] } else { rc[f + 1] } * u;
            }

            // 3) monta as linhas (laço sem dependência entre iterações → vetoriza)
            for i in 1..n - 1 {
                let (ge, gd) = (g[i - 1], g[i]);
                let (fe, fd) = (fl[i - 1], fl[i]);
                let cap = rc[i] * cap_fator;
                aw[i] = -ge - fe.max(0.0);
                ap[i] = cap + ge + gd + fd.max(0.0) + (-fe).max(0.0) - (fd - fe);
                ae[i] = -gd - (-fd).max(0.0);
                bb[i] = cap * t[i];
            }

            // 4) Thomas "twisted": elimina de cima pra baixo E de baixo pra cima no mesmo laço.
            //    São duas cadeias de dependência independentes → a CPU executa as duas em
            //    paralelo e a latência da divisão é escondida pela metade.
            //    Os nós 0 e n-1 entram como linhas identidade (= condições de contorno).
            //    Cima (linhas 1..=m):   x_i + c_i x_{i+1} = d_i
            //    Baixo (linhas m+1..):  x_j + e_j x_{j-1} = f_j   (guardados nos mesmos vetores)
            let (mut c, mut d) = (0.0, t[0]);
            let (mut e, mut f) = (0.0, t[n - 1]);
            for s in 0..h {
                let i = 1 + s;
                let inv = 1.0 / (ap[i] - aw[i] * c);
                c = ae[i] * inv;
                d = (bb[i] - aw[i] * d) * inv;
                c_l[i] = c;
                d_l[i] = d;

                let j = n - 2 - s;
                let inv = 1.0 / (ap[j] - ae[j] * e);
                e = aw[j] * inv;
                f = (bb[j] - ae[j] * f) * inv;
                c_l[j] = e;
                d_l[j] = f;
            }
            if m > h {
                // número ímpar de linhas: a metade de cima tem uma a mais
                let inv = 1.0 / (ap[m] - aw[m] * c);
                c_l[m] = ae[m] * inv;
                d_l[m] = (bb[m] - aw[m] * d) * inv;
            }

            // 5) encontro no meio: resolve o 2x2 das linhas m e m+1
            let (cm, dm, em, fm) = (c_l[m], d_l[m], c_l[m + 1], d_l[m + 1]);
            let mut xc = (dm - cm * fm) / (1.0 - cm * em);
            let mut xb = fm - em * xc;
            t[m] = xc;
            t[m + 1] = xb;

            // 6) volta pros dois lados, também intercalada
            for s in 1..h {
                let i = m - s;
                xc = d_l[i] - c_l[i] * xc;
                t[i] = xc;
                let j = m + 1 + s;
                xb = d_l[j] - c_l[j] * xb;
                t[j] = xb;
            }
            if m > h {
                t[1] = d_l[1] - c_l[1] * xc;
            }
        }

        self.meio1.temperaturas.copy_from_slice(&t[..n1]);
        self.meio2.temperaturas.copy_from_slice(&t[n1..]);
    }

    fn plotar(&self, arquivo: &str, titulo: &str) {
        let root = BitMapBackend::new(arquivo, (800, 600)).into_drawing_area();

        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .caption(titulo, ("sans-serif", 30))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..self.metodo.dimensao, 0.0..350.0)
            .unwrap();

        chart
            .configure_mesh()
            .x_desc("Posição [m]")
            .y_desc("Temperatura [°C]")
            .draw()
            .unwrap();

        let temperaturas = self
            .meio1
            .temperaturas
            .iter()
            .chain(self.meio2.temperaturas.iter());

        chart
            .draw_series(LineSeries::new(
                temperaturas.enumerate().map(|(i, &temp)| {
                    let x = i as f64 * self.metodo.dx;
                    (x, temp)
                }),
                &RED,
            ))
            .unwrap();

        root.present().unwrap();
    }
}

//  Em Rust, arrays têm tamanho fixo na compilação. Para um tamanho definido em
// tempo de execução, utiliza-se  Vec<f32> :

//Se o tamanho for fixo e conhecido previamente, use const generics:

fn criar_array<const N: usize>(valor: f64) -> [f64; N] {
    [valor; N]
}

// a pira aqui é passar o meio, a posiçao do "t1" e um Option caso seja a interface. também é
// importante passar os parâmetros do método, no caso , dx e a area onde passa o fluxo que vai ser 1
// sempre eu acho
// parametros antigos:ta: f64, tb: f64, ka: f64, kb: f64, dxa: f64, dxb: f64, area: f64
fn interface_flux(sistema: &Sistema, posicao: usize) -> (f64, f64) {
    let (meio1, indice1) = sistema.indice(posicao);
    let (meio2, indice2) = sistema.indice(posicao + 1);

    let dx = sistema.metodo.dx;

    let t1 = meio1.temperaturas[indice1];
    let t2 = meio2.temperaturas[indice2];

    let k1 = meio1.material.get_k(t1);
    let k2 = meio2.material.get_k(t2);

    let r = dx / (k1 * 2.0) + dx / (k2 * 2.0); //resistencia
    let g = 1.0 / r; //condutancia termica
    let q = g * (t1 - t2); // fluxo 
    (q, g)
}

//esboço

// aqui vai a função que junta os dois métodos e compara eles lado a lado de acordo com o
//número de iterações no tempo
//
//entrada: sistema implicito e sistema explicito. também o número de tempo
fn comparacao(sisimpl: &Sistema, sisexpl: &Sistema) {
    //
    //tirar as diferenças

    //printar os dados relevantes.
    for i in 0..sisimpl.meio1.temperaturas.len() {
        if i % 100 == 0 {
            println!(
                "oia a dif: {:?} na posiçao  {:?}",
                sisexpl.meio1.temperaturas[i] - sisimpl.meio1.temperaturas[i],
                i
            )
        }
    }
    for i in 0..sisimpl.meio2.temperaturas.len() {
        if i % 100 == 0 {
            println!(
                "oia a dif: {:?} na posiçao  {:?}",
                sisexpl.meio2.temperaturas[i] - sisimpl.meio2.temperaturas[i],
                i + sisexpl.meio1.temperaturas.len()
            )
        }
    }
}

fn maior_diferenca(sisimpl: &Sistema, sisexpl: &Sistema) {
    let mut maior: f64 = 0.0;
    let mut posicao: usize = 0;
    for i in 0..sisimpl.meio1.temperaturas.len() {
        let dif = (sisexpl.meio1.temperaturas[i] - sisimpl.meio1.temperaturas[i]).abs();
        if dif > maior {
            maior = dif;
            posicao = i;
        }
    }
    println!("maior dif no meio1: {:?} posição: {:?}", maior, posicao);
    for i in 0..sisimpl.meio2.temperaturas.len() {
        let dif = (sisexpl.meio2.temperaturas[i] - sisimpl.meio2.temperaturas[i]).abs();
        if dif > maior {
            maior = dif;
            posicao = i + sisexpl.meio1.temperaturas.len();
        }
    }
    println!("maior dif global: {:?} posição: {:?}", maior, posicao);
}

fn benchmark<F: FnOnce()>(iteracoes: usize, funcao: F) {
    let inicio = Instant::now();

    funcao();

    let tempo = inicio.elapsed().as_secs_f64();

    println!("Tempo de execução: {:.4} s", tempo);

    if iteracoes > 0 {
        println!(
            "Tempo por iteração: {:.4} µs",
            tempo * 1e6 / iteracoes as f64
        );
    }
}

#[inline]
fn preencher_props(material: &Material, temps: &[f64], k: &mut [f64], rc: &mut [f64]) {
    match material {
        Material::Ar(tabela) => {
            for ((&t, k), rc) in temps.iter().zip(k.iter_mut()).zip(rc.iter_mut()) {
                let p = tabela.consultar(t);
                *k = p.k;
                *rc = p.rho * p.cp;
            }
        }
        m => {
            for ((&t, k), rc) in temps.iter().zip(k.iter_mut()).zip(rc.iter_mut()) {
                *k = m.get_k(t);
                *rc = m.get_rho(t) * m.get_cp(t);
            }
        }
    }
}
