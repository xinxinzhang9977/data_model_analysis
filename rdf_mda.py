import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
import matplotlib.pyplot as plt


def calculate_average_rdf(topology, trajectory, selection_A, selection_B,
                          total_frames=100, frame_interval=10,
                          r_min=2,r_max=15.0, nbins=100):
    """
    计算平均RDF的函数

    参数:
    topology: 拓扑文件路径
    trajectory: 轨迹文件路径
    selection_A, selection_B: 原子选择语句
    total_frames: 总帧数
    frame_interval: 计算间隔帧数
    r_max: 最大径向距离
    nbins: 分段数
    """

    # 加载轨迹
    try:
        universe = mda.Universe(topology, trajectory)

        print(f"成功加载轨迹: {trajectory}")
    except Exception as e:
        print(f"加载轨迹失败: {e}")
        return None, None

    # 选择原子组
    try:
        group_A = universe.select_atoms(selection_A)
        group_B = universe.select_atoms(selection_B)
        print(f"组A原子数: {len(group_A)}")
        print(f"组B原子数: {len(group_B)}")
    except Exception as e:
        print(f"原子选择失败: {e}")
        return None, None

    # 检查总帧数
    actual_frames = len(universe.trajectory)
    print(f"轨迹实际帧数: {actual_frames}")

    if total_frames > actual_frames:
        total_frames = actual_frames
        print(f"调整总帧数为: {total_frames}")

    # 存储所有RDF结果
    all_rdf_results = []
    rdf_bin_centers = None

    # 分段计算RDF
    for start_frame in range(0, total_frames, frame_interval):
        end_frame = min(start_frame + frame_interval - 1, total_frames - 1)

        print(f"计算帧 {start_frame}-{end_frame} 的RDF...")

        # 创建RDF分析对象
        rdf = InterRDF(group_A, group_B,
                       nbins=nbins,
                       range=(r_min, r_max))

        # 计算当前区间的RDF
        rdf.run(start=start_frame, stop=end_frame, step=1)

        # 存储结果
        all_rdf_results.append(rdf.results.rdf)

        # 保存bin中心位置（只需一次）
        if rdf_bin_centers is None:
            rdf_bin_centers = rdf.results.bins

    # 计算平均RDF
    average_rdf = np.mean(all_rdf_results, axis=0)

    return rdf_bin_centers, average_rdf


# 使用示例
if __name__ == "__main__":
    # 文件路径 - 请根据你的实际情况修改
    filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\5.MD_e_6V/'
    topology_file = filepath+'NVEe.data'  # 如: 'system.pdb', 'system.gro'等
    trajectory_file = filepath+'NVEe_1_5.lammpsdump'  # 如: 'md.xtc', 'md.dcd'等

    # 原子选择 - 请根据你的实际情况修改
    selection_mg = 'resid 801 to 820'  # 例如选择Mg原子
    selection_dme = 'resid 1 to 800'  # 例如选择DME
    selection_tfsi = 'resid 821 to 860' # 选择TFSI
    nbins = 500
    total_frames = 1000
    frame_interval = 100
    # 计算平均RDF
    bins, avg_rdf_dme = calculate_average_rdf(
        topology=topology_file,
        trajectory=trajectory_file,
        selection_A=selection_mg,
        selection_B=selection_dme,
        total_frames=total_frames,
        frame_interval=frame_interval,
        r_min=2.0,
        r_max=15.0,
        nbins=nbins
    )
    bins, avg_rdf_tfsi = calculate_average_rdf(
        topology=topology_file,
        trajectory=trajectory_file,
        selection_A=selection_mg,
        selection_B=selection_tfsi,
        total_frames=total_frames,
        frame_interval=frame_interval,
        r_min=2.0,
        r_max=15.0,
        nbins=nbins
    )
    bins, avg_rdf_mg = calculate_average_rdf(
        topology=topology_file,
        trajectory=trajectory_file,
        selection_A=selection_mg,
        selection_B=selection_mg,
        total_frames=total_frames,
        frame_interval=frame_interval,
        r_min=2.0,
        r_max=15.0,
        nbins=nbins
    )

    if bins is not None and avg_rdf_dme is not None and avg_rdf_tfsi is not None:
        # 设置全局字体大小
        plt.rcParams.update({
            'font.size': 18,
            'axes.titlesize': 18,
            'axes.labelsize':18,
            'xtick.labelsize': 18,
            'ytick.labelsize': 18,
            'legend.fontsize': 18,
            'figure.dpi': 100
        })

        # 创建图形
        fig, ax = plt.subplots(figsize=(8,6))

        # 绘制曲线
        ax.plot(bins, avg_rdf_dme, 'r-', linewidth=2, label='DME')
        ax.plot(bins, avg_rdf_tfsi, 'b-', linewidth=2, label='TFSI')
        ax.plot(bins, avg_rdf_mg, 'g-', linewidth=2, label='Mg')

        # 设置坐标轴标签和标题
        ax.set_xlabel('Radial Distance (Å)', fontweight='bold')
        ax.set_ylabel('g(r)', fontweight='bold')
        ax.set_title('Average Radial Distribution Function', fontweight='bold', pad=15)

        # 设置图例（去掉边框）
        ax.legend(frameon=False)

        # 设置网格
        ax.grid(alpha=0.3)

        # 设置坐标轴边框加粗
        for spine in ax.spines.values():
            spine.set_linewidth(2)

        # 设置刻度朝内
        ax.tick_params(axis='both',
                       direction='in',  # 刻度朝内
                       length=6,  # 刻度线长度
                       width=2,  # 刻度线宽度
                       top=True,  # 显示上侧刻度
                       right=True,  # 显示右侧刻度
                       labelsize=18)  # 刻度标签大小

        # 添加次要刻度
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1)

        # 设置刻度数字加粗
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('normal')

        # 调整布局
        plt.tight_layout()
        plt.show()

        # 保存数据
        output_name = trajectory_file + '_average_rdf_results.dat'
        np.savetxt(output_name,
                   np.column_stack((bins, avg_rdf_dme, avg_rdf_tfsi, avg_rdf_mg)),
                   header='Distance(A)\tg_dme(r)\tg_tfsi(r)\tg_mg(r)',
                   fmt='%.6f',
                   delimiter='\t')
        print(f"结果已保存到 {output_name}")