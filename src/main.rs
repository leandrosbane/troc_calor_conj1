use plotters::prelude::*;

fn main() {
    println!("Hello, world!");

    // definindo variáveis
    // passar tudo issopra enum
    //let ferrinho = "material_a";
    //let cobrezinho = "material_b";

    //definindo a malha de pontos
    let tamanho: f64 = 1.0;
    const N: usize = 1000;
    const ITERACOES_TEMPO: usize = 1_000_000;
    let n: f64 = N as f64;
    let dx = tamanho / n;
    let dt = 0.01;

    //condições de contorno
    let ta = 0.0;
    let tb = 350.;

    //meios
    let mut meio1 = Meio::criar_material(N / 2, 0.0, Material::Ferro);
    let mut meio2 = Meio::criar_material(N / 2, 0.0, Material::Cobre);

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
    let mut sistema = Sistema {
        meio1,
        meio2,
        metodo: metodinho,
        interface: N / 2 - 1,
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
    sistema.calc_explicito(ITERACOES_TEMPO);

    sistema.plotar("explicito.png", "Método explícito");

    //parte do implicito
    //let coeficientes = sistema.coef_implicit(6);
    //println!("Coeficientes do nó 6: {:?}", coeficientes);
    //sistema.met_implicit();

    //loop completo
    let mut tridiag = SistemaTridiagonal::criar(N - 2);

    println!("len b = {}", tridiag.b.len());
    println!("destino meio1 = {}", sistema.meio1.temperaturas[1..].len());

    println!("{:?}", &tridiag.diagonal_principal[..5]);
    println!("{:?}", &tridiag.b[..5]);

    for _ in 0..ITERACOES_TEMPO {
        sistema.met_implicit(&mut tridiag);

        sistema.solver_implicit(&mut tridiag);

        // atualizar meio1 e meio2 com temperaturas_internas
        sistema.meio1.temperaturas[1..].copy_from_slice(&tridiag.b[..499]);

        sistema.meio2.temperaturas[..499].copy_from_slice(&tridiag.b[499..]);
    }

    //vetor completo
    let mut temperaturas = Vec::with_capacity(N);

    temperaturas.extend_from_slice(&sistema.meio1.temperaturas);
    temperaturas.extend_from_slice(&sistema.meio2.temperaturas);

    //plot
    sistema.plotar("implicito.png", "Método implícito");
}

struct Meio {
    temperaturas: Vec<f64>,
    material: Material,
}

struct Metodo {
    dimensao: f64,
    nos: usize,
    dx: f64,
    dt: f64,
}

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

enum Material {
    Ferro,
    Cobre,
    Brazino,
}

impl Material {
    fn get_k(&self, t: f64) -> f64 {
        match self {
            Material::Ferro => 10.0 + 0.5 * t,
            Material::Cobre => 15.0 + 0.3 * t,
            Material::Brazino => 20.0 + 0.1 * t,
        }
    }

    fn get_rho(&self) -> f64 {
        match self {
            Material::Ferro => 7874.0,
            Material::Cobre => 8960.0,
            Material::Brazino => 5000.0,
        }
    }

    fn get_cp(&self) -> f64 {
        match self {
            Material::Ferro => 449.0,
            Material::Cobre => 385.0,
            Material::Brazino => 500.0,
        }
    }
}

impl Meio {
    fn criar_material(tamanho: usize, valor: f64, material: Material) -> Self {
        Self {
            temperaturas: vec![valor; tamanho],
            material: material,
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

            let rho = mei.material.get_rho();
            let cp = mei.material.get_cp();
            let termo_sei_la_o_nome = self.metodo.dt / (rho * cp * self.metodo.dx);

            let atualtemp = mei.temperaturas[pos];

            let novatemp = termo_sei_la_o_nome * (q_esquerda - q_direita) + atualtemp;

            q_esquerda = q_direita;

            if j < q1 {
                buffer1[pos] = novatemp
            } else {
                buffer2[pos] = novatemp
            }
        }
        //TODO acho que a parada hoje é mudar a saída do código pra dois vetores separados.
    }

    fn coef_implicit(&self, indice_global: usize) -> (f64, f64, f64, f64) {
        //essa função só vai retornar os coeficientes pra usar em cada posição da matriz
        let (meio, indice_local) = self.indice(indice_global);

        let rho = meio.material.get_rho();
        let cp = meio.material.get_cp();
        let capacidade = (rho * cp * self.metodo.dx) / self.metodo.dt;

        let (_, ge) = interface_flux(self, indice_global - 1);
        let (_, gd) = interface_flux(self, indice_global);

        let temp_atual = meio.temperaturas[indice_local];

        //aW, aP, aE, b
        (-ge, capacidade + ge + gd, -gd, capacidade * temp_atual)
    }

    fn calc_explicito(&mut self, num_iteracoes: usize) {
        let mut novas_temperaturas_meio1 = self.meio1.temperaturas.clone();
        let mut novas_temperaturas_meio2 = self.meio2.temperaturas.clone();

        for i in 0..num_iteracoes {
            self.iteracao(&mut novas_temperaturas_meio1, &mut novas_temperaturas_meio2);

            std::mem::swap(&mut self.meio1.temperaturas, &mut novas_temperaturas_meio1);
            std::mem::swap(&mut self.meio2.temperaturas, &mut novas_temperaturas_meio2);

            //colocar um if nesse println pra não fazer 1 MELHAO de prints!!
            let passo = i + 1;
            if passo % 100000 == 0 {
                println!(
                    "iteração {:?}, começo: {:?} fim: {:?}",
                    i,
                    &self.meio1.temperaturas[0..5],
                    &self.meio2.temperaturas[495..500],
                );
                //print da interface
                println!("Fe: {:?}", &self.meio1.temperaturas[495..500]);
                println!("Cu: {:?}", &self.meio2.temperaturas[0..5]);
            }
        }
    }

    //função funcionando
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
            ) = self.coef_implicit(i + 1);
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

    //thompsoooooooooonnnnnnnnnnnnnnn
    fn solver_implicit(&self, sistema: &mut SistemaTridiagonal) {
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
            sistema.b[n - 1] = (sistema.b[n - 1] - sistema.diagonal_superior[n - 1] * sistema.b[n])
                / sistema.diagonal_principal[n - 1];
            sistema.diagonal_principal[n - 1] = 1.0;

            //repete até o indice zero
        }
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
